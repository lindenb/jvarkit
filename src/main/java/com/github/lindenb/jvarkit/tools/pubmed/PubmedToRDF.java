/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.Maps;


/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="pubmed2rdf",
	description="Programming language use distribution from recent programs / articles",
	keywords={"pubmed","xml","rdf","gene"},
	creationDate="20250528",
	modificationDate="20250528",
	generate_doc = false,
	jvarkit_amalgamion =  true,
	menu="Pubmed"
	)
public class PubmedToRDF
	extends Launcher
	{
	private static final String EVI="https://w3id.org/EVI#";
	private static final String DCTERMS="http://purl.org/dc/terms/";
	private static final Logger LOG = Logger.of(PubmedToRDF.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;
	@Parameter(names={"--rec"},description="rec file")
	private List<Path> recFiles=new ArrayList<>();
	
	
	private final Map<String,EWord> entities_nodes = new HashMap<>(50_000);
	private final Map<String,String> gene2uri = new HashMap<>(50_000);
	
	private static class Known {
		Set<String> uri;
	}
	
	private static abstract class ENode {
		EWord parent;
		ENode firstChild;
		ENode nextSibling;
		abstract void scan(List<String> words,int idx,Consumer<String> consumer);
		String getPath() {
			List<String> L=new ArrayList<String>();
			EWord p = this.parent;
			while(p!=null) {
				L.add(p.text);
				p=p.parent;
				}
			Collections.reverse(L);
			return String.join(" ", L);
			}
		}
	
	private static class EWord extends ENode {
		final String text;
		EWord(final String text) {
			this.text = text;
			}
		EWord find(final String s) {
			ENode n= this.firstChild;
			while(n!=null)
				{
				if((n instanceof EWord)) {
					EWord w = EWord.class.cast(n);
					if(w.text.equalsIgnoreCase(s)) return w;
					}
				n=n.nextSibling;
				}
			return null;
			}
		
		ENode append(ENode n) {
			n.parent=this;
			if(this.firstChild==null) {
				firstChild = n;
				}
			else
				{
				n.nextSibling=this.firstChild;
				this.firstChild=n;
				}
			return n;
			}
		void add(final List<String> words,int idx,String uri) {
			if(idx==words.size()) {
				final EURI euri = new EURI();
				euri.uri = uri;
				append(euri);
				}
			else {
				EWord w = find(words.get(idx)); 
				if(w==null) {
					w=new EWord(words.get(idx));
					append(w);
					}
					
				w.add(words,idx+1,uri);
				}
			}

		
		@Override
		void scan(List<String> words,int idx,Consumer<String> consumer) {
			if(idx >= words.size()) return;
			boolean b = words.get(idx).equalsIgnoreCase(this.text);
			if(!b) return;
			ENode c = this.firstChild;
			while(c!=null) {
				c.scan(words, idx+1, consumer);
				c=c.nextSibling;
				}
			}
		
		@Override
		public String toString() {
			StringBuilder sb=new StringBuilder(this.text);
			sb.append("{");
			ENode c = this.firstChild;
			while(c!=null) {
				sb.append(c.toString());
				c=c.nextSibling;
				}
			sb.append("}");
			return sb.toString();
			}
		}
	
	private static class EURI extends ENode {
		String uri;
		void scan(List<String> words,int idx,Consumer<String> consumer) {
			consumer.accept(uri);
			LOG.info("GOT uri = "+uri+" "+getPath());
			}
		@Override
		public String toString() {
			return "{URI="+uri+"}";
			}
		}
	
	
	
	private static class Entry implements Comparable<Entry> {
		String pmid;
		String title="";
		String pubDate;
		final Set<String> entities=new HashSet<>();
		@Override
		public int compareTo(Entry o) {
			int i= this.pubDate.compareTo(o.pubDate);
			if(i!=0) return i;
			return this.pmid.compareTo(o.pmid);
		}
		
		void write(XMLStreamWriter w) throws XMLStreamException {
			w.writeStartElement("evi","Article",EVI);
			w.writeAttribute("rdf", RDF.NS, "about", "https://pubmed.ncbi.nlm.nih.gov/"+pmid);
			
			w.writeStartElement("dcterms","title",DCTERMS);
			w.writeCharacters(this.title);
			w.writeEndElement();
			
			w.writeStartElement("dcterms","issued",DCTERMS);
			w.writeCharacters(this.pubDate);
			w.writeEndElement();

			w.writeEndElement();
		}
	}
		
		
	public PubmedToRDF()
		{
		
		}
	
	private void skip(final XMLEventReader r) throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				skip(r);
				}
			}
		}
	
	
	/** extract text from stream. Cannot use XMLEventReader.getTextContent() 
	 * when a public title contains some tag like '<sup>'
	 */
	private String textContent(final XMLEventReader r) throws XMLStreamException
		{
		final StringBuilder sb=new StringBuilder();
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				sb.append(textContent(r));
				}
			else if(evt.isCharacters())
				{
				sb.append(evt.asCharacters().getData());
				}
			}
		return sb.toString();
		}
	
	private List<String> tokenizeSentences(final String text) {
		return Arrays.asList(text.split("\\. "));
	}
	
	
	private String ScanDate(final XMLEventReader r) throws XMLStreamException {
		String year = null;
		String month="01";
		String day="01";
		final Map<String,String> month2month= Maps.of(
				"Jan","01","Feb","02","Mar","03","Apr","04","May","05","Jun","06",
				"Jul","07","Aug","08","Sep","09","Oct","10","Nov","11","Dec","12"
				);
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement())
				{
				String name =evt.asStartElement().getName().getLocalPart();
				if(name.equals("Year")) {
					year = r.getElementText();
					}
				else if(name.equals("Month")) {
					month =   r.getElementText();
					month = month2month.getOrDefault(month, month);
					}
				else if(name.equals("Day")) {
					day = r.getElementText();
					}
				else
					{
					skip(r);
					}
				}
			else if(evt.isEndElement()) {
				if(year!=null) {
					return String.join("-", year,month,day);
					}
				return null;
				}
			}
		
		 return null;
		}
	
	private Entry scanArticle(
			final String rootName,
			final XMLEventReader r
			) throws XMLStreamException,IOException {
			final Entry entry=new Entry();
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				
				if(evt.isStartElement()) {
					final StartElement start = evt.asStartElement();
					final String eltName = start.getName().getLocalPart();
					if(entry.pmid==null && eltName.equals("PMID")) {
						entry.pmid=r.getElementText();
					}
					else if(entry.title.isEmpty() && eltName.equals("ArticleTitle")) {
						entry.title= textContent(r);
					}
					else if(eltName.equals("PubDate")) {
						entry.pubDate = ScanDate(r);
					}
					else if(eltName.equals("NameOfSubstance")) {
						String subst= r.getElementText().trim();
						final String protein_human=", protein, human";
						int idx=subst.indexOf(protein_human);
						if(idx!=-1 && (idx+protein_human.length())==subst.length()) {
							final String gene_name = subst.substring(0,idx);
							String uri = this.gene2uri.get(gene_name);
							if(uri!=null) entry.entities.add(uri);
							}
						}
					else if(eltName.equals("Abstract")) {
						for(String abstractText : tokenizeSentences(textContent(r))) {
							List<String> words = Arrays.asList(abstractText.split("[ ]+"));
							for(int i=0;i< words.size();i++) {
								final String gene_name = words.get(i);
								final String uri = this.gene2uri.get(gene_name);
								if(uri!=null) {
									entry.entities.add(uri);
									}
								
								EWord eword= this.entities_nodes.get(words.get(i).toLowerCase());
								if(eword!=null) {
									eword.scan(words, i, S->entry.entities.add(S));
									}
								}
							}
						}
					}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals(rootName)) 
					{
					if(!entry.entities.isEmpty()) {
						
						return entry;
						}
					return null;
					}
				}
			
			}//end of xml read
		return null;
		}
	
	private void readRecFile(final Path path) throws IOException {
		try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
			final List<Map.Entry<String, String>> items=new ArrayList<>();
			for(;;) {
				String line= br.readLine();
				if(StringUtils.isBlank(line)) {
					if(!items.isEmpty()) {
						String type= items.stream().
								filter(KV->KV.getKey().equals("type")).
								map(KV->KV.getValue()).
								collect(HtsCollectors.toSingleton());
						
						if(type.equals("gene")) {
							String symbol= items.stream().filter(KV->KV.getKey().equals("symbol")).
									map(KV->KV.getValue()).
									collect(HtsCollectors.toSingleton());
							String uri= items.stream().filter(KV->KV.getKey().equals("uri")).
									map(KV->KV.getValue()).
									collect(HtsCollectors.toSingleton());
							this.gene2uri.put(symbol, uri);
							}
						else if(type.equals("disease")) {
							
							String uri= items.stream().filter(KV->KV.getKey().equals("uri")).
									map(KV->KV.getValue()).
									collect(HtsCollectors.toSingleton());
							for(String text:  items.stream().filter(KV->KV.getKey().equals("text")).
										map(KV->KV.getValue()).
										collect(Collectors.toSet())) {
								final List<String> words = Arrays.asList(text.replaceAll("[ ., ]+"," ").split("[ ]+"));
								final String startw= words.get(0).toLowerCase();
								EWord eword = this.entities_nodes.get(startw);
								if(eword==null) {
									eword = new EWord(startw);
									this.entities_nodes.put(startw, eword);
									}
								eword.add(words,1,uri);
								LOG.info("adding "+eword);
								}
							}
						else
							{
							throw new IllegalArgumentException("unknown type "+type);
							}
						}
					if(line==null) break;
					items.clear();
					continue;
					}
				int colon=line.indexOf(':');
				if(colon==-1) throw new IOException("colon missing in "+line);
				String k = line.substring(0,colon).trim().toLowerCase();
				String v = line.substring(colon+1).trim();
				if(!StringUtils.isBlank(k) && !StringUtils.isBlank(v)) {
					items.add(new AbstractMap.SimpleEntry<>(k, v));
					}
			}
		}
	}
	
	
	@Override
	public int doWork(final List<String> args) {
	
		try {
			final XMLInputFactory xmlInputFactory = XMLInputFactory.newFactory();
			xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.debug("Ignoring resolve Entity");
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			final String inputName= oneFileOrNull(args);
			for(Path recPath: this.recFiles) {
				readRecFile(recPath);
				}
			final List<Entry> entries=new ArrayList<>();
			try(InputStream in=(inputName==null?stdin():IOUtils.openURIForReading(inputName))) {
				final XMLEventReader r = xmlInputFactory.createXMLEventReader(in);
				
				while(r.hasNext()) {
					final XMLEvent evt= r.nextEvent();
					if(evt.isStartElement() )	{
						final String localName= evt.asStartElement().getName().getLocalPart();
						if(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle"))
							{
							Entry entry=scanArticle(localName,r);
							if(entry!=null) {
								entries.add(entry);
								Collections.sort(entries);
								while(entries.size()>100) {
									entries.remove(entries.size()-1);
									}
								}
							}
						}
					}
				
				try(PrintStream out = super.openPathOrStdoutAsPrintStream(this.outFile)) {
					final XMLOutputFactory xof = XMLOutputFactory.newFactory();
					final XMLStreamWriter w = xof.createXMLStreamWriter(out, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					w.writeStartElement("rdf:RDF");
					w.writeNamespace("rdf", RDF.NS);
					w.writeNamespace("rdf", RDF.NS);
					
					for(Entry entry: entries) {
						entry.write(w);
						}
					
					w.writeEndElement(); //
					w.writeEndDocument();
					w.close();
					out.flush();
					}
				
				r.close();
				}
			
			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
		} 

		}
		
	
	public static void main(String[] args) {
		new PubmedToRDF().instanceMainWithExit(args);
	}
}
