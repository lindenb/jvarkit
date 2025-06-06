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
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ncbi.PubmedUtils;
import com.github.lindenb.jvarkit.rdf.RDFModel;
import com.github.lindenb.jvarkit.rdf.Resource;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;
import com.github.lindenb.jvarkit.stream.HtsCollectors;


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
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private static final Resource GENE_RSRC= new Resource("http://purl.obolibrary.org/obo/","SO_0000704");
	private static final Resource DISEASE_RSRC = new Resource("http://purl.obolibrary.org/obo/","DOID_4");
	private static final Resource ANNOT_PROP = new Resource(NS,"annotation");
	private static final Resource ARTICLE_TYPE = new Resource(NS,"Article");
	private static final Resource PUB_DATE_PROP = new Resource("http://purl.org/dc/terms/","date");
	private static final Resource TITLE_PROP = new Resource("http://purl.org/dc/terms/","title");
	private static final Logger LOG = Logger.of(PubmedToRDF.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;
	@Parameter(names={"--rec"},description="rec file")
	private List<Path> recFiles=new ArrayList<>();
	
	
	private final Map<String,EWord> entities_nodes = new HashMap<>(50_000);
	private final Map<String,Resource> gene2hgnc = new HashMap<>(50_000);
	
	private static class Known {
		Set<String> uri;
		}
	
	private static abstract class ENode {
		EWord parent;
		ENode firstChild;
		ENode nextSibling;
		abstract void scan(List<String> words,int idx,final Entry entry);
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
		
		void add(final List<String> words,int idx,final Resource uri, final String disease_name) {
			if(idx==words.size()) {
				final EDisease euri = new EDisease(uri,disease_name);
				append(euri);
				}
			else {
				EWord w = find(words.get(idx)); 
				if(w==null) {
					w=new EWord(words.get(idx));
					append(w);
					}
					
				w.add(words,idx+1,uri,disease_name);
				}
			}

		
		@Override
		void scan(List<String> words,int idx,final Entry entry) {
			if(idx >= words.size()) return;
			boolean b = words.get(idx).equalsIgnoreCase(this.text);
			if(!b) return;
			ENode c = this.firstChild;
			while(c!=null) {
				c.scan(words, idx+1, entry);
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
	
	private static class EDisease extends ENode {
		private final Resource uri;
		private final String name;
		public EDisease(final Resource uri, final String name) {
			this.uri = uri;
			this.name = name;
			}
		
		void scan(final List<String> words,int idx,final Entry entry) {
			entry.addDisease(name, uri);
			}
		@Override
		public String toString() {
			return "{URI="+uri+"}";
			}
		}
	
	
	
	private static class Entry implements Comparable<Entry> {
		Resource pmid;
		final RDFModel statements = new RDFModel();
		
		public Date getPubDate() {
			return this.statements.
					findMatching(pmid, PUB_DATE_PROP, null).
					map(it->it.getObject().asLiteral().getDate()).
					findAny().
					orElse(null);
			}
		
		@Override
		public int compareTo(final Entry o) {
			final int i= this.getPubDate().compareTo(o.getPubDate());
			if(i!=0) return i;
			return this.pmid.getLocalName().compareTo(o.pmid.getLocalName());
			}
		
		void addGene(String geneName,Resource hgnc) {
			this.statements.addStatement(
				hgnc,
				RDF.type,
				GENE_RSRC
				);
			this.statements.addStatement(
					hgnc,
				RDFS.label,
				geneName
				);
			
			this.statements.addStatement(
				this.pmid,
				ANNOT_PROP,
				hgnc
				);
			}
		
		void addDisease(String diseaseName,Resource uri) {
			this.statements.addStatement(
				uri,
				RDF.type,
				DISEASE_RSRC
				);
			this.statements.addStatement(
				uri,
				RDFS.label,
				diseaseName
				);
			this.statements.addStatement(
				this.pmid,
				ANNOT_PROP,
				uri
				);
			}
		
		
	}
		
		
	public PubmedToRDF()
		{
		
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
						final String pmid =r.getElementText();
						entry.pmid=new Resource(PubmedUtils.pmidToURL(pmid));
						entry.statements.addStatement(entry.pmid, RDF.type, ARTICLE_TYPE);
						}
					else if(entry.pmid!=null && entry.statements.findMatching(entry.pmid,TITLE_PROP,null).count()==0L && eltName.equals("ArticleTitle")) {
						String title= textContent(r);
						entry.statements.addStatement(entry.pmid, TITLE_PROP, title);
					}
					else if(entry.pmid!=null && eltName.equals("PubDate")) {
						final Optional<Date> pubDate = PubmedUtils.scanDate(r);
						if(pubDate.isPresent()) entry.statements.addStatement(entry.pmid, PUB_DATE_PROP, pubDate.get());
					}
					else if(eltName.equals("NameOfSubstance")) {
						String subst= r.getElementText().trim();
						final String protein_human=", protein, human";
						int idx=subst.indexOf(protein_human);
						if(idx!=-1 && (idx+protein_human.length())==subst.length()) {
							final String gene_name = subst.substring(0,idx);
							final Resource hgnc = this.gene2hgnc.get(gene_name);
							if(hgnc!=null) {
								entry.addGene(gene_name, hgnc);
								}
							}
						}
					else if(eltName.equals("Abstract")) {
						for(String abstractText : tokenizeSentences(textContent(r))) {
							List<String> words = Arrays.asList(abstractText.split("[ ]+"));
							for(int i=0;i< words.size();i++) {
								final String gene_name = words.get(i);
								final Resource hgnc = this.gene2hgnc.get(gene_name);
								if(hgnc!=null) {
									entry.addGene(gene_name,hgnc);
									}
								
								final EWord eword= this.entities_nodes.get(words.get(i).toLowerCase());
								if(eword!=null) {
									eword.scan(words, i, entry);
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
					if(!entry.statements.isEmpty()) {
						
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
							final Resource hgnc= items.stream().filter(KV->KV.getKey().equals("hgnc")).
									map(KV->KV.getValue()).
									map(S->new Resource("http://identifiers.org/hgnc/",S)).
									collect(HtsCollectors.toSingleton());
							this.gene2hgnc.put(symbol, hgnc);
							}
						
						else if(type.equals("disease")) {
							Resource uri= items.stream().filter(KV->KV.getKey().equals("uri")).
									map(KV->KV.getValue()).
									map(S->new Resource(S)).
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
								eword.add(words,1,uri,text);
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
				
				final RDFModel model=new RDFModel();
				entries.stream().forEach(E->model.addAll(E.statements));
								
				try(PrintStream out = super.openPathOrStdoutAsPrintStream(this.outFile)) {
					model.writeXml(out);
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
