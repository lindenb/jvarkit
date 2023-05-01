/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.psi;

import java.util.ArrayList;
import java.util.List;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

public class AbstractPsiObject {
protected static final QName QNAME_ID = new QName("id");
protected static class Names  {
	String shortLabel;
	String fullName;
	}

protected static class Reference  {
	boolean primary;
	String db;
	String id;
	}
	
protected String id = null;
protected Names names = null;
protected final List<Reference> references = new ArrayList<>();

public String getId() {
	return this.id;
	}
@Override
public int hashCode() {
	return getId().hashCode();
	}

public Names getNames() {
	return this.names;
}


protected static Names parseNames(XMLEventReader r) throws XMLStreamException {
	final Names n = new Names();
	while(r.hasNext()) {
		final XMLEvent evt = r.nextEvent();
		if(evt.isStartElement()) {
			final String lcl  = evt.asStartElement().getName().getLocalPart();
			if(lcl.equals("shortLabel")) {
				n.shortLabel = r.getElementText();
				}
			else if(lcl.equals("fullName")) {
				n.fullName = r.getElementText();
				}
			}
		else if(evt.isEndElement()) {
			final String lcl  = evt.asEndElement().getName().getLocalPart();
			if(lcl.equals("names")) break;
			}
		}
	return n;
	}
protected static List<Reference> parseXref(XMLEventReader r) throws XMLStreamException {
	final List<Reference> L = new ArrayList<>();
	while(r.hasNext()) {
		final XMLEvent evt = r.nextEvent();
		if(evt.isStartElement()) {
			final StartElement startE = evt.asStartElement();
			final String lcl  = startE.getName().getLocalPart();
			if(lcl.equals("primaryRef") || lcl.equals("secondaryRef")) {
				Attribute att = startE.getAttributeByName(new QName("db"));
				Reference ref = new Reference();
				ref.primary = lcl.charAt(0)=='p';
				ref.db = att.getValue();
				att = startE.getAttributeByName(QNAME_ID);
				ref.id = att.getValue();
				L.add(ref);
				}
			}
		else if(evt.isEndElement()) {
			final String lcl  = evt.asEndElement().getName().getLocalPart();
			if(lcl.equals("xref")) break;
			}
		}
	return L;
	}
}
