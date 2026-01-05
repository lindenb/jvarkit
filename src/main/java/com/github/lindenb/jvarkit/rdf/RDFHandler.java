/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
*/
package com.github.lindenb.jvarkit.rdf;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.rdf.ns.RDF;



/**
 * RDFHandler
 * RDF reader calling the method 'found' each time it finds a
 * new RDFStatement
 *
 */
public class RDFHandler {
	private static final Logger LOG = Logger.of(RDFHandler.class);
	private static final QName RDF_about = new QName(RDF.NS,"about");
	private final XMLInputFactory factory;
	private Resource base=null;
	private final Consumer<Statement> consumer;

	
	/**
	 * constructor
	 * initialize the XMLInputFactory
	 */
	public RDFHandler(final Consumer<Statement> consumer)
		{
		this.factory = XMLInputFactory.newInstance();
		this.factory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
		this.factory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
		this.factory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
		this.consumer = consumer;
		}
	

	private void found(Statement stmt) {
		this.consumer.accept(stmt);
	}
	
	
	
	public void setBase(Resource base) {
		this.base = base;
		}

	

	
	/** do the parsing */
	private void read(Reader in) throws XMLStreamException
		{
		try {
			/** create a XML Event parser */
			XMLEventReader parser = this.factory.createXMLEventReader(in);
			/** loop until we find a rdf:RDF element */
			while(parser.hasNext())
				{
				XMLEvent event = parser.nextEvent();
				if(event.isStartElement())
					{
					StartElement start=(StartElement)event;
					LOG.debug("found start "+start);
					if(name2string(start).equals(RDF.NS+"RDF"))
						{
						LOG.debug("found RDF");
						parseRDF(parser);
						}
					}
				}
			}
		catch (URISyntaxException e)
			{
			throw new XMLStreamException(e);
			}
		}

	
/** we're in the rdf:RDF , scan all statements */
private void parseRDF(final XMLEventReader reader) throws XMLStreamException,URISyntaxException
	{
	LOG.debug("In RDF ROOT: loop over each rdf:Description");
	while(reader.hasNext())
		{
		final XMLEvent event = reader.nextEvent();
		if(event.isEndElement())
			{
			//we're done
			return;
			}
		else if(event.isStartElement())
			{
			parseDescription(reader,event.asStartElement());
			}
		else
			{
			otherEvent(event);
			}
		}
	}

private void otherEvent(XMLEvent event) throws XMLStreamException {
	 if(event.isProcessingInstruction())
		{
		throw new XMLStreamException("Found Processing Instruction in RDF ???",event.getLocation());
		}
	else if(event.isCharacters() &&
			!StringUtils.isBlank(event.asCharacters().getData()))
		{
		throw new XMLStreamException("Found text in RDF ???",event.getLocation());
		}
	else
		{
		throw new XMLStreamException("boum",event.getLocation());
		}
	}

/**
 * Parse description of a Resource
 * @param description
 * @return
 * @throws URISyntaxException
 * @throws XMLStreamException
 */
private Resource parseDescription(final XMLEventReader reader,final StartElement description) throws URISyntaxException,XMLStreamException
	{
	LOG.debug("Found a new  rdf:Description "+description.getName());
	final Resource subject;
	Attribute att= description.getAttributeByName(RDF_about);
	if(att!=null) {
		subject= createURI( att.getValue());
		}
	else
		{
		att= description.getAttributeByName(new QName(RDF.NS,"nodeID"));
		if(att!=null) {
			subject=  createURI( att.getValue());
			}
		else
			{
			att= description.getAttributeByName(new QName(RDF.NS,"ID"));
			if(att!=null) {
				subject= resolveBase(att.getValue());
				}
			else
				{
				subject= createAnonymousURI();
				}
			}
		}
	
	
	LOG.debug("Description uri=\""+subject+"\"");
	
	QName qn=description.getName();
	if(!(qn.getNamespaceURI().equals(RDF.NS) &&
		 qn.getLocalPart().equals("Description")))
		{
		found(new Statement(subject,RDF.type, new Resource(qn)));
		}

	/** loop over attributes */
	for(Iterator<?> i=description.getAttributes();
		i.hasNext();)
		{
		att=(Attribute)i.next();
		qn= att.getName();
		final String local= qn.getLocalPart();
		if(qn.getNamespaceURI().equals(RDF.NS) &&
			( local.equals("about") ||
				local.equals("ID") ||
				local.equals("nodeID")))
				{
				continue;
				}
		
		found(new Statement(subject,new Resource(qn),new Literal(att.getValue())));
		}
	
	
	while(reader.hasNext())
		{
		final XMLEvent event = reader.nextEvent();
		
		if(event.isEndElement())
			{
			return subject;
			}
		else if(event.isStartElement())
			{
			parsePredicate(reader,subject,event.asStartElement());
			}
		else 
			{
			otherEvent(event);
			}
		}

	return subject;
	}

/**
 * parse predicate
 * @param descriptionURI
 * @param predicate 
 * @throws URISyntaxException
 * @throws XMLStreamException
 */
private void parsePredicate(XMLEventReader reader,Resource descriptionURI,StartElement predicate) throws URISyntaxException,XMLStreamException
	{
	String parseType=null;
	String lang=null;
	Resource datatype=null;
	Attribute att;
	QName qn=null;
	Resource resource=null;

	Resource predicateURI=new Resource(predicate.getName());
	LOG.debug("parse rdf:description=\""+descriptionURI+"\" predicate:"+predicateURI);
	
	/** collect attributes */
	for(int loop=0;loop<2;++loop)
		{
		for(Iterator<?> i=predicate.getAttributes();
			i.hasNext();)
			{
			att=(Attribute)i.next();
			qn= att.getName();
			String local= qn.getLocalPart();
			if(qn.getPrefix().equals("xml") &&
				local.equals("lang"))
				{
				if(loop==0) lang=att.getValue();
				continue;
				}
			else if(qn.getNamespaceURI().equals(RDF.NS))
				{
				if(local.equals("parseType"))
					{
					if(loop==0)  parseType=att.getValue();
					LOG.debug("parseType:"+parseType);
					continue;
					}
				else if(local.equals("datatype"))
					{
					if(loop==0) datatype= createURI(att.getValue());
					LOG.debug("dataType="+datatype);
					continue;
					}
				else if(local.equals("resource"))
					{
					if(loop==0) resource=  createURI(att.getValue());
					LOG.debug("rdf:resource="+resource);
					continue;
					}
				else if(local.equals("nodeID"))
					{
					if(loop==0) resource=  createURI(att.getValue());
					LOG.debug("rdf:nodeID="+resource);
					continue;
					}
				else if(local.equals("ID"))
					{
					if(loop==0) resource= resolveBase(att.getValue());
					LOG.debug("ID="+resource);
					continue;
					}
				}
			
			if(loop==1)
				{
				if(resource!=null)
					{
					found(new Statement(
							resource,
							new Resource(att.getName()),
							new Literal(att.getValue())
							));
					}
				else
					{
					throw new XMLStreamException("Cannot handle attribute "+att);
					}
				}
			
			}
		}
	
	
	
	
	if(resource!=null)
		{
		found(new Statement(descriptionURI, predicateURI, resource));
		XMLEvent event= reader.peek();
		if(event!=null && event.isEndElement())
			{
			reader.nextEvent();
			return;
			}
		throw new XMLStreamException("Expected a EndElement for this element");
		}
	
	if(parseType==null) parseType="default";
	
	if(parseType.equals("Literal"))
		{
		final StringBuilder b= parseLiteral(reader);
		found(new Statement(
				descriptionURI,
				predicateURI,
				parseLiteral(b.toString(), lang, datatype))
				);
		}
	else if(parseType.equals("Resource"))
		{
		final Resource blanck = createAnonymousURI();
		found(new Statement(descriptionURI, predicateURI, blanck));
		while(reader.hasNext())
			{
			final XMLEvent event = reader.nextEvent();
			if(event.isStartElement())
				{
				parsePredicate(reader,blanck, event.asStartElement());
				}
			else if(event.isEndElement())
				{
				return;
				}
			else
				{
				otherEvent(event);
				}
			}
		
		}
	else  if(parseType.equals("Collection"))
		{
		final List<Statement> items=new ArrayList<>();
		while(reader.hasNext())
			{
			final XMLEvent event = reader.nextEvent();
			if(event.isStartElement())
				{
				Resource value= parseDescription(reader,event.asStartElement());
				items.add(new Statement(
						descriptionURI,
						predicateURI,
						value
						));
				
				}
			else if(event.isEndElement())
				{
				for(Statement i: items) {
					// todo convert to list
					found(i);
					}
				return;
				}
			else
				{
				otherEvent(event);
				}
			}
		}
	else
		{
		boolean foundResourceAsChild=false;
		StringBuilder b= new StringBuilder();
		while(reader.hasNext())
			{
			XMLEvent event = reader.nextEvent();
			if(event.isStartElement())
				{
				if(b.toString().trim().length()!=0)
					{
					throw new XMLStreamException(
							"Bad text \""+b+"\" before "+
							event.asStartElement().getName()
							);
					}
				final Resource childURI=parseDescription(reader,event.asStartElement());
				found(new Statement(descriptionURI,predicateURI,childURI));
				b.setLength(0);
				foundResourceAsChild=true;
				}
			else if(event.isCharacters())
				{
				b.append(event.asCharacters().getData());
				}
			else if(event.isEndElement())
				{
				if(!foundResourceAsChild)
					{
					found(new Statement(
						descriptionURI,
						predicateURI,
						parseLiteral(b.toString(),lang,datatype)
						));
					}
				else
					{
					if(b.toString().trim().length()!=0) throw new XMLStreamException("Found bad text "+b);
					}
				return;
				}
			}
		}
	
}

private Resource resolveBase(String ID) throws URISyntaxException
	{
	if(this.base==null) return  createURI(ID);
	return new Resource(this.base.toURI().resolve(ID));
	}

private Literal parseLiteral(String text,String lang,Resource xsdType) {
	return new Literal(text);
	}

private StringBuilder parseLiteral(XMLEventReader reader) throws XMLStreamException
	{
	StringBuilder b=new StringBuilder();
	QName qn;
	int depth=0;
	while(reader.hasNext())
		{
		XMLEvent event = reader.nextEvent();
		if(event.isCharacters())
			{
			b.append(escapeXML(event.asCharacters().getData()));
			}
		else if(event.isProcessingInstruction())
			{
			b.append("<?"+event.asCharacters()+"?>");
			}
		else if(event.isEndElement())
			{
			if(depth==0) return b;
			qn= event.asEndElement().getName();
			b.append("</"+qn.getPrefix()+":"+qn.getLocalPart()+">");
			depth--;
			}
		else if(event.isStartElement())
			{
			qn= event.asEndElement().getName();
			b.append("<"+qn.getPrefix()+":"+qn.getLocalPart());
			
			for(Iterator<?> i=event.asStartElement().getAttributes();
				i.hasNext();)
				{
				Attribute att=(Attribute)i.next();
				qn= att.getName();
				b.append(" ").
					append(qn.getPrefix()+":"+qn.getLocalPart()).
					append("=\"").
					append(escapeXML(att.getValue())).
					append("\"");
				}
			event= reader.peek();
			if(event!=null && event.isEndElement())
				{
				reader.nextEvent();
				b.append("/>");
				}
			else
				{
				b.append(">");
				depth++;
				}
			}
		}
	
	return b;
	}

protected Resource createAnonymousURI() throws URISyntaxException
	{
	return new Resource();
	}

public void parse(InputStream in) throws XMLStreamException
	{
	read(new InputStreamReader(in));
	}

public void parse(Reader in) throws XMLStreamException
	{
	read(in);
	}

public void parse(final Path in) throws XMLStreamException
	{
	try(Reader fin=IOUtils.openFileForReader(in)) {
		read(fin);
		}
	catch (IOException e)
		{
		throw new XMLStreamException(e);
		}
	}

public void parse(URL in) throws XMLStreamException
	{
	LOG.debug("parsing URL "+in);
	try {
		try(InputStream fin= in.openStream()) {
			read(new InputStreamReader(fin));
			}
		}
	catch (IOException e)
		{
		throw new XMLStreamException(e);
		}
	}

private String name2string(StartElement e)
	{
	return name2string(e.getName());
	}

private String name2string(QName name)
	{
	return name.getNamespaceURI()+name.getLocalPart();
	}


private Resource createURI(String uriAsString) throws URISyntaxException 
	{
	return  new Resource(uriAsString);
	}


private static String escapeXML(String s)
	{
	return StringUtils.escapeXML(s);
	}

}
