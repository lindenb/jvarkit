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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.file.Path;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.rdf.ns.OWL;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;

import htsjdk.samtools.util.DateParser;

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
	private static final Logger LOG = Logger.of(PubmedToRDF.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;
	@Parameter(names={"--mondo"},description="Mondo ontology URI")
	private String mondoURI="https://purl.obolibrary.org/obo/mondo.owl";
	
	
	private final Set<Entity> entities = new HashSet<>();
	
	private static class Entry implements Comparable<Entry> {
		String pmid;
		String title="";
		Calendar pubDate;
		final Set<Entity> entities=new HashSet<>();
		@Override
		public int compareTo(Entry o) {
			int i= this.pubDate.compareTo(o.pubDate);
			if(i!=0) return i;
			return this.pmid.compareTo(o.pmid);
		}
		
		void write(XMLStreamWriter w) throws XMLStreamException {
			w.writeStartElement("entry");
			
			w.writeStartElement("title");
			w.writeCharacters(this.title);
			w.writeEndElement();
			
			w.writeEmptyElement("link");
			w.writeAttribute("href", ""+pmid);

			
			w.writeStartElement("id");
			w.writeCharacters(UUID.nameUUIDFromBytes(pmid.getBytes()).toString());
			w.writeEndElement();
			
			
			w.writeStartElement("content");
			w.writeAttribute("type", "text/html");
			
			w.writeEndElement();
			
			w.writeEndElement();
		}
	}
	
	private abstract class Entity {
		final String uri;
		Entity(final String uri) {
			this.uri = uri;
			}
		public String getURI() {
			return uri;
			}
		abstract boolean find(List<String> words);
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null || !(obj instanceof Entity))  return false;
			return uri.equals(Entity.class.cast(obj).uri);
			}
		@Override
		public int hashCode() {
			return uri.hashCode();
			}
		}
	
	private class GeneEntity extends Entity {
			final String symbol;
			
			GeneEntity(final String symbol,final String uri) {
				super(uri);
				this.symbol = symbol;
				}
			
			@Override
			boolean find(List<String> words) {
				return words.stream().anyMatch(S->S.equals(this.symbol));
				}
			
			@Override
			public String toString() {
				return symbol;
				}
			
		 }
	
	private class DiseaseEntity extends Entity {
		final List<String> tokens;
		DiseaseEntity(final String label,final String uri) {
			super(uri);
			this.tokens = null;//TODO
			}
		@Override
		boolean find(List<String> words) {
			if(tokens.size()> words.size()) return false;
			for(int i=0;i+tokens.size()<= words.size();i++) {
				int x=0;
				while(x< this.tokens.size()) {
					if(!tokens.get(x).equalsIgnoreCase(words.get(i+x))) {
						break;
						}
					}
				if(x==tokens.size()) {
					return true;
					}
				}
			return false;
			}

		@Override
		public String toString() {
			return String.join(" ", tokens);
			}
		
		}
	
	
	
	
	
		
	public PubmedToRDF()
		{
		
		}
	private static final QName OWL_Class = new QName(OWL.NS,"Class");
	private static final QName RDFS_label = new QName(RDFS.NS,"label");
	private static final QName RDF_about = new QName(RDF.NS,"about");
	
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
	

	private  void visitOWLClass(XMLEventReader r,StartElement root) throws XMLStreamException {
		final Attribute att=root.getAttributeByName(RDF_about);
		if(att==null) {
			skip(r);
			return;
			}
		String uri = att.getValue();
		String label=null;
		while(r.hasNext()) {
			XMLEvent evt=r.nextTag();
			if(evt.isStartElement()) {
				final StartElement se  = evt.asStartElement();
				if(label==null && se.getName().equals(RDFS_label)) {
					label = r.getElementText();
					}
				else
					{
					skip(r);
					}
				}
			else if(evt.isEndElement()) {
				if(uri.startsWith("http://identifiers.org/hgnc/") && !StringUtils.isBlank(label)) {
					final GeneEntity entity=new GeneEntity(label,uri);
					entities.add(entity);
					}
				else
					{
					final  DiseaseEntity entity = new DiseaseEntity(label, uri);
					entities.add(entity);
					}
				
				return;
				}
			}
		}
	
	private void scanMondoOntoly(String uri) throws IOException,XMLStreamException {
		XMLInputFactory  xmlInputFactory = XMLInputFactory.newFactory();
		xmlInputFactory.setXMLResolver(new XMLResolver() {
			@Override
			public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
					throws XMLStreamException {
				LOG.debug("Ignoring resolve Entity");
				return new ByteArrayInputStream(new byte[0]);
			}
		});
		try(InputStream in= IOUtils.openURIForReading(uri)) {
			final XMLEventReader r=xmlInputFactory.createXMLEventReader(in);
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				if(evt.isStartElement()) {
					final StartElement E=evt.asStartElement();
					if(E.getName().equals(OWL_Class)) {
						visitOWLClass(r, E);
						}
					}
				}
			r.close();
			}
		LOG.info("number of entities :"+this.entities.size());
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
		return Arrays.asList(text.split("(?<=[.!?])\\s+"));
	}
	
	private List<String> tokenizeWords(final String text) {
		return null;
		}
	
	private Calendar ScanDate(final XMLEventReader r) throws XMLStreamException {
		String year = null;
		String month="01";
		String day="01";
		
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
			String article_year=null;
			boolean PubDate=false;
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
					else if(article_year==null && PubDate && eltName.equals("Year")) {
						article_year=r.getElementText();
					}
					else if(eltName.equals("NameOfSubstance")) {
						String subst= r.getElementText().trim();
						final String protein_human=", protein, human";
						int idx=subst.indexOf(protein_human);
						if(idx!=-1 && (idx+protein_human.length())==subst.length()) {
							final String gene_name = subst.substring(0,idx);
							this.entities.stream().
								filter(G->(G instanceof GeneEntity)).
								map(G->GeneEntity.class.cast(G)).
								filter(G->G.symbol.equals(gene_name)).
								findFirst().
								ifPresent(G->entry.entities.add(G));
							}
						}
					else if(eltName.equals("Abstract")) {
						for(String abstractText : tokenizeSentences(textContent(r))) {
							List<String> words = tokenizeWords(abstractText);
							for(Entity entity:this.entities) {
								if(entity.find(words)) {
									entry.entities.add(entity);
									}
								}
							}
						}
					}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals("PubDate")) {
					PubDate=false;
					}
				else if(eltName.equals(rootName)) 
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
			this.scanMondoOntoly(this.mondoURI);
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
					w.writeStartElement("feed");
					w.writeDefaultNamespace("http://www.w3.org/2005/Atom");
					
					for(Entry entry: entries) {
						entry.write(w);
						}
					
					w.writeEndElement();
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
