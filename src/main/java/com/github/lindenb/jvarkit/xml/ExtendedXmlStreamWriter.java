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

*/
package com.github.lindenb.jvarkit.xml;

import java.io.Closeable;
import java.io.Flushable;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Objects;

import javax.xml.namespace.NamespaceContext;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.w3c.dom.Comment;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.RuntimeIOException;

public class ExtendedXmlStreamWriter implements Closeable, Flushable {
	private final XMLStreamWriter delegate;

	public ExtendedXmlStreamWriter(final XMLStreamWriter delegate) {
		this.delegate = Objects.requireNonNull(delegate);
		}
	
	public XMLStreamWriter getDelegate() {
		return delegate;
		}
	
	private void error(final XMLStreamException err) {
		throw new RuntimeIOException(err);
		}
	
	private String toString(Object o) {
		return String.valueOf(o);
	}
	
	public ExtendedXmlStreamWriter writeStartElement(QName qName)  {
		if(StringUtils.isBlank(qName.getNamespaceURI())) {
			return writeStartElement(qName.getLocalPart());
			}
		
		if(StringUtils.isBlank(qName.getPrefix())) {
			return writeStartElement(qName.getNamespaceURI(),qName.getLocalPart());
			}
		return writeStartElement(qName.getPrefix(), qName.getLocalPart(), qName.getNamespaceURI());
		}
	
	public ExtendedXmlStreamWriter writeElement(QName qName,Object value)  {
		writeStartElement(qName);
		writeCharacters(value);
		return writeEndElement();
		}
	
	public ExtendedXmlStreamWriter writeEmptyElement(QName qName)  {
		return writeEmptyElement(qName.getPrefix(), qName.getLocalPart(), qName.getNamespaceURI());
		}
	
	
	public ExtendedXmlStreamWriter writeElement(String localName,Object value)  {
		writeStartElement(localName);
		writeCharacters(value);
		return writeEndElement();
		}
	
	public AutoCloseable element(final String localName) {
		writeStartElement(localName);
		return new AutoCloseable() {
			@Override
			public void close() throws Exception {
				writeEndElement();
				}
			};
		}
	
	public ExtendedXmlStreamWriter writeStartElement(String localName)  {
		try {
			getDelegate().writeStartElement(localName);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeElement(String namespaceURI, String localName,Object value) {
		writeStartElement(namespaceURI, localName);
		writeCharacters(value);
		writeEndElement();
		return this;
		}
	
	
	public ExtendedXmlStreamWriter writeStartElement(String namespaceURI, String localName) {
		try {
			getDelegate().writeStartElement(namespaceURI,localName);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public AutoCloseable element(String namespaceURI, String localName) {
		writeStartElement(namespaceURI,localName);
		return new AutoCloseable() {
			@Override
			public void close() throws Exception {
				writeEndElement();
				}
			};
		}

	
	
	public ExtendedXmlStreamWriter writeElement(String prefix, String localName, String namespaceURI,Object value) {
		writeStartElement(prefix, localName, namespaceURI);
		writeCharacters(value);
		return writeEndElement();
		}

	
	public ExtendedXmlStreamWriter writeStartElement(String prefix, String localName, String namespaceURI)  {
		try {
			getDelegate().writeStartElement(prefix, localName, namespaceURI);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEmptyElement(String namespaceURI, String localName) {
		try {
			getDelegate().writeEmptyElement(namespaceURI, localName);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEmptyElement(String prefix, String localName, String namespaceURI) {
		try {
			getDelegate().writeEmptyElement(prefix, localName, namespaceURI);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEmptyElement(String localName) {
		try {
			getDelegate().writeEmptyElement(localName);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEndElement() {
		try {
			getDelegate().writeEndElement();
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEndDocument() {
		try {
			getDelegate().writeEndDocument();
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	@Override
	public void close() {
		try {
			getDelegate().close();
		} catch(XMLStreamException err) {
			error(err);
			}
		}

	@Override
	public void flush() {
		try {
			getDelegate().flush();
		} catch(XMLStreamException err) {
			error(err);
			}
		}

	public ExtendedXmlStreamWriter writeAttributes(final Map<String,Object> atts) {
		for(final String key:atts.keySet()) {
			writeAttribute(key, atts.get(key));
			}
		return this;
		}
	
	public ExtendedXmlStreamWriter writeAttribute(String localName, Object value) {
		try {
			getDelegate().writeAttribute(localName,toString(value));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeAttribute(String prefix, String namespaceURI, String localName, Object value) {
		try {
			getDelegate().writeAttribute(prefix, namespaceURI, localName,toString(value));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeAttribute(String namespaceURI, String localName, Object value) {
		try {
			getDelegate().writeAttribute(namespaceURI, localName,toString(value));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}
	
	public ExtendedXmlStreamWriter writeAttribute(final QName qName, Object value) {
		return writeAttribute(qName.getPrefix(), qName.getNamespaceURI(), qName.getLocalPart(), value);
		}

	public ExtendedXmlStreamWriter writeNamespace(String prefix, String namespaceURI) {
		try {
			getDelegate().writeNamespace(prefix, namespaceURI);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeDefaultNamespace(String namespaceURI) {
		try {
			getDelegate().writeDefaultNamespace(namespaceURI);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeComment(Object data) {
		try {
			getDelegate().writeComment(toString(data));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeProcessingInstruction(String target) {
		try {
			getDelegate().writeProcessingInstruction(target);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeProcessingInstruction(String target, Object data)  {
		try {
			getDelegate().writeProcessingInstruction(target,toString(data));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeCData(Object data)  {
		try {
			getDelegate().writeCData(toString(data));
		} catch(XMLStreamException err) {
			error(err);
			}
		return  this;
		}

	public ExtendedXmlStreamWriter writeDTD(String dtd)  {
		try {
			getDelegate().writeDTD(dtd);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeEntityRef(String name)  {
		try {
			getDelegate().writeEntityRef(name);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeStartDocument()  {
		try {
			getDelegate().writeStartDocument();
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeStartDocument(String version)   {
		try {
			getDelegate().writeStartDocument(version);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeStartDocument(String encoding, String version) {
		try {
			getDelegate().writeStartDocument(encoding, version);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeStartDocument(Charset encoding, String version) {
		return writeStartDocument(encoding.displayName(), version);
		}
	
	
	public ExtendedXmlStreamWriter nl() {
		return writeCharacters("\n");
	}
	
	
	public ExtendedXmlStreamWriter writeDOM(Node root) {
		if(root==null) return this;
		switch(root.getNodeType()) {
			case Node.DOCUMENT_NODE:
				writeStartDocument(Document.class.cast(root).getXmlVersion(),Document.class.cast(root).getInputEncoding());
				for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling()) writeDOM(n1);
				return writeEndDocument();
			case Node.DOCUMENT_FRAGMENT_NODE:
				for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling()) writeDOM(n1);
				return this;
			case Node.COMMENT_NODE:
				return writeComment(Comment.class.cast(root).getData());
			case Node.TEXT_NODE:
				return writeCharacters(org.w3c.dom.Text.class.cast(root).getData());
			case Node.CDATA_SECTION_NODE:
				return writeCData(org.w3c.dom.CDATASection.class.cast(root).getData());
			case Node.ATTRIBUTE_NODE:
				return writeAttribute( root.getPrefix(),root.getNamespaceURI(), root.getLocalName(),root.getNodeValue());
			case Node.ELEMENT_NODE:
				if(!root.hasChildNodes()) {
					writeEmptyElement(root.getPrefix(), root.getLocalName(), root.getNamespaceURI());
					if(root.hasAttributes()) {
						final NamedNodeMap nm = root.getAttributes();
						for(int i=0;i< nm.getLength();i++) writeDOM(nm.item(i));
						}
					}
				else
					{
					writeStartElement(root.getPrefix(), root.getLocalName(), root.getNamespaceURI());
					if(root.hasAttributes()) {
						final NamedNodeMap nm = root.getAttributes();
						for(int i=0;i< nm.getLength();i++) writeDOM(nm.item(i));
						}
					for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling()) writeDOM(n1);
					writeEndDocument();
					}
				return this;
			default: throw new IllegalStateException();
			}
		}
	
	public ExtendedXmlStreamWriter writeCharacters(Object o) {
		try {
			getDelegate().writeCharacters(toString(o));
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter writeCharacters(char[] text, int start, int len)  {
		try {
			getDelegate().writeCharacters(text,start,len);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public String getPrefix(String uri)  {
		try {
			return getDelegate().getPrefix(uri);
		} catch(XMLStreamException err) {
			error(err);
			throw new IllegalStateException(err);
			}
		}

	public ExtendedXmlStreamWriter setPrefix(String prefix, String uri)  {
		try {
			getDelegate().setPrefix(prefix, uri);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter setDefaultNamespace(String uri) {
		try {
			getDelegate().setDefaultNamespace(uri);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public ExtendedXmlStreamWriter setNamespaceContext(NamespaceContext context)  {
		try {
			getDelegate().setNamespaceContext(context);
		} catch(XMLStreamException err) {
			error(err);
			}
		return this;
		}

	public NamespaceContext getNamespaceContext() {
		return getDelegate().getNamespaceContext();
	}

	public Object getProperty(String name) {
		return getDelegate().getProperty(name);
	}

}
