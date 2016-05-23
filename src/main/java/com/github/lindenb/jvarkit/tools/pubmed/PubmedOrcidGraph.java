/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;


import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

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
import com.google.gson.stream.JsonWriter;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.CloserUtil;


public class PubmedOrcidGraph
	extends AbstractPubmedOrcidGraph
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(PubmedOrcidGraph.class);

	
	private static class Author implements Comparable<Author>
		{
		String foreName = null;
		String lastName = null ;
		String orcid = null;
		String initials=null;
		String affiliation =null;
		
		/** idea: in xml+xslt , used to presort orcid indentifiers
		 * in order to produce unique pairs of collab(orcid1,orcid2) */
		@Override
		public int compareTo(final Author o) {
			return this.orcid.compareTo(o.orcid);
			}
		
		void json(final JsonWriter w) throws IOException {
			w.beginObject();
			w.name("orcid");w.value(this.orcid);
			if(foreName!=null) {
				w.name("foreName");w.value(this.foreName);
			}
			if(lastName!=null) {
				w.name("lastName");w.value(this.lastName);
			}
			if(initials!=null) {
				w.name("initials");w.value(this.initials);
			}
			if(affiliation!=null) {
				w.name("affiliation");w.value(this.affiliation);
			}
			w.endObject();
		}
		
		void xml(final XMLStreamWriter w) throws XMLStreamException {
			w.writeStartElement("Author");
			w.writeAttribute("orcid", orcid);
			
			if(foreName!=null) {
				w.writeStartElement("foreName");
				w.writeCharacters(foreName);
				w.writeEndElement();
			}
			
			if(lastName!=null) {
				w.writeStartElement("lastName");
				w.writeCharacters(lastName);
				w.writeEndElement();
			}
			
			if(initials!=null) {
				w.writeStartElement("initials");
				w.writeCharacters(initials);
				w.writeEndElement();
			}
			if(affiliation!=null) {
				w.writeStartElement("affiliation");
				w.writeCharacters(affiliation);
				w.writeEndElement();
			}
			w.writeEndElement();
		}
		
		void print(PrintWriter pw) {
			pw.println("#Author");
			pw.println("ORCID\t"+orcid);
			if(foreName!=null) pw.println("ForeName\t"+foreName);
			if(lastName!=null) pw.println("LastName\t"+lastName);
			if(initials!=null) pw.println("Initials\t"+initials);
			if(affiliation!=null) pw.println("Affiliation\t"+affiliation);
			pw.println();
		}
		}
	
	private Author parseAuthor(XMLEventReader r,final String pmid)  throws XMLStreamException
		{
		final Author au = new Author();
		while(r.hasNext()) {
			final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				String eltName = start.getName().getLocalPart();
				if(eltName.equals("LastName")) {
					au.lastName=r.getElementText().trim();
				} else if(eltName.equals("ForeName")) {
					au.foreName=r.getElementText().trim();
				} else if(eltName.equals("Affiliation")) {
					au.affiliation=r.getElementText().trim();
				} else if(eltName.equals("Initials")) {
					au.initials=r.getElementText().trim();
				} else if(eltName.equals("Identifier")) {
					final Attribute source= start.getAttributeByName(new QName("Source"));
					
					if(source!=null && !source.getValue().equalsIgnoreCase("ORCID")) {
						LOG.warn("interesting, a non-orcid Identifier in pmid:"+pmid+" "+source.getValue());
						}
					else if(source!=null && source.getValue().equalsIgnoreCase("ORCID")) {
						au.orcid=r.getElementText().trim();
						final int slash = au.orcid.lastIndexOf('/');
						if(slash!=-1) au.orcid=au.orcid.substring(slash+1);
						au.orcid=au.orcid.trim();
						
						if(au.orcid.startsWith("orcid.org.")) {
							LOG.warn("Ugly orcid in "+au.orcid +" pmid:"+ pmid);
							final int dot= au.orcid.lastIndexOf('.');
							au.orcid=au.orcid.substring(dot+1);
						}
						else if(au.orcid.startsWith("orcid.org") ) {
							LOG.warn("Ugly orcid in "+au.orcid +" pmid:"+ pmid);
							au.orcid=au.orcid.substring(9);
						}
						
						
						if(au.orcid.length()==16){
							// 0000000226481522 instead of 0000-0002-2648-1522
							au.orcid= au.orcid.substring(0, 4 ) +"-"+
									  au.orcid.substring(4, 8 ) +"-"+
									  au.orcid.substring(8, 12) +"-"+
									  au.orcid.substring(12)
									  ;
						}
					}
				}
			}
			else if(evt.isEndElement() &&
				evt.asEndElement().getName().getLocalPart().equals("Author")) {
				return au;
			}
		}
		throw new IllegalStateException("should never happen");
		}
	
	
	
	private void scanArticle(
			final String rootName,
			final XMLEventReader r,
			final PrintWriter pw,
			final XMLStreamWriter w,
			final JsonWriter jsw
			) throws XMLStreamException,IOException {
		try {
			String pmid=null;
			String ArticleTitle=null;
			String Year=null;
			String ISOAbbreviation=null;
			String doi  = null;
			boolean PubDate=false;
			List<Author> authors = new ArrayList<>();
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				final String eltName = start.getName().getLocalPart();
				if(pmid==null && eltName.equals("PMID")) {
					pmid=r.getElementText();
				}
				else if(ISOAbbreviation==null && eltName.equals("ISOAbbreviation")) {
					ISOAbbreviation=r.getElementText();
				}
				else if(ArticleTitle==null && eltName.equals("ArticleTitle")) {
					ArticleTitle=r.getElementText();
				}
				else if(eltName.equals("PubDate")) {
					PubDate=true;
				}
				else if(doi==null && eltName.equals("ArticleId")) {
					final Attribute idType= start.getAttributeByName(new QName("IdType"));
					if(idType!=null && idType.getValue().equalsIgnoreCase("doi")) {
						doi = r.getElementText().trim();
					}

				}
				else if(Year==null && PubDate && eltName.equals("Year")) {
					Year=r.getElementText();
				}
				else if(eltName.equals("Author")) {
					final Author author = parseAuthor(r,pmid);
					if(author.orcid!=null) authors.add(author);
				}
				}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals("PubDate")) {
					PubDate=false;
				}
				if(eltName.equals(rootName)) break;
			}
			}//end of xml read
		Collections.sort(authors);
		if(authors.isEmpty()) {
			//do nothing
		}
		else if(w!=null) {
			w.writeStartElement(rootName);
			w.writeAttribute("pmid", pmid);
			if(doi!=null) w.writeAttribute("doi", doi);
			if(Year!=null)
			{
				w.writeStartElement("year");
				w.writeCharacters(Year);
				w.writeEndElement();
			}
			if(ISOAbbreviation!=null)
			{
				w.writeStartElement("journal");
				w.writeCharacters(ISOAbbreviation);
				w.writeEndElement();
			}
			if(ArticleTitle!=null)
			{
				w.writeStartElement("title");
				w.writeCharacters(ArticleTitle);
				w.writeEndElement();
			}
			for(final Author au:authors) {
				au.xml(w);
			}
			w.writeEndElement();
		}
		else if(jsw!=null) {
			jsw.beginObject();
			
			jsw.name("pmid");
			jsw.value(pmid);
			if(doi!=null){
				jsw.name("doi");
				jsw.value(doi);
			}
			if(Year!=null)
			{
				jsw.name("year");
				jsw.value(Year);
			}
			if(ISOAbbreviation!=null)
			{
				jsw.name("journal");
				jsw.value(ISOAbbreviation);
			}
			if(ArticleTitle!=null)
			{
				jsw.name("title");
				jsw.value(ArticleTitle);
			}
			
			jsw.name("authors");
			jsw.beginArray();
				for(final Author au:authors) {
					au.json(jsw);
				}
			jsw.endArray();
			jsw.endObject();
		}
		else
			{
			pw.print("#");
			pw.println(rootName);
			pw.println("PMID\t"+pmid);
			if(doi!=null) pw.println("doi\t"+doi);
			if(Year!=null) pw.println("Year\t"+Year);
			if(ISOAbbreviation!=null) pw.println("Journal\t"+ISOAbbreviation);
			if(ArticleTitle!=null) pw.println("Title\t"+ArticleTitle);
			for(final Author au:authors) {
				pw.println("Author\t"+au.orcid);
			}
			pw.println();
			
			for(final Author au:authors) {
				au.print(pw);
			}
			}
			
			
		
		} finally {
		}
	}
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		InputStream in=null;
		XMLEventReader r = null;
		XMLStreamWriter w=null;
		PrintWriter pw=null;
		JsonWriter jsw=null;
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
			in=(inputName==null?stdin():IOUtils.openURIForReading(inputName));
			r = xmlInputFactory.createXMLEventReader(in);
			
			pw= super.openFileOrStdoutAsPrintWriter();
			if("xml".equalsIgnoreCase(super.format)) {
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				w = xof.createXMLStreamWriter(pw);
				w.writeStartDocument( "UTF-8","1.0");
				w.writeStartElement("PubmedArticleSet");
				w.writeComment("Generated with "+getName()+" "+getOnlineDocUrl()+" - Pierre Lindenbaum.");
			} else if("json".equalsIgnoreCase(super.format)) {
				jsw = new JsonWriter(pw);
				jsw.beginArray();
			}
			while(r.hasNext()) {
				final XMLEvent evt= r.nextEvent();
				if(evt.isStartElement() )	{
					final String localName= evt.asStartElement().getName().getLocalPart();
					if(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle"))
						{
						scanArticle(localName,r,pw,w,jsw);
						}
				}
			}
			
			if(jsw!=null) {
				jsw.endArray();
				jsw.close();
			}
			if(w!=null) {
				w.writeEndElement();
				w.writeEndDocument();
				w.flush();w.close();
				w=null;
				}
			pw.flush();
			pw.close();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(r);
			CloserUtil.close(pw);
			CloserUtil.close(w);
		}
		}

	
	public static void main(String[] args)
		{
		new PubmedOrcidGraph().instanceMainWithExit(args);
		}
	}
