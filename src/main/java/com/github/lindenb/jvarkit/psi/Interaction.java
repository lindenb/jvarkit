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


public class Interaction extends AbstractPsiObject {

private final List<Participant> participants = new ArrayList<>();
	
private Interaction() {
}
	

static class Participant {
	String interactorId=null;
	String role = null;
	}

@Override
public boolean equals(final Object obj) {
	if(obj==this) return true;
	if(obj==null || !(obj instanceof Interaction)) return false;
	return this.getId().equals(Interaction.class.cast(obj).getId());
	}

private static Participant parseParticipant(XMLEventReader r) throws XMLStreamException {
	final Participant p = new Participant();
	while(r.hasNext()) {
		final XMLEvent evt = r.nextEvent();
		if(evt.isStartElement()) {
			final StartElement startE = evt.asStartElement();
			final String lcl  = startE.getName().getLocalPart();
			if(lcl.equals("proteinInteractorRef")) {
				final Attribute att = startE.getAttributeByName(new QName("ref"));
				if(att!=null) p.interactorId = att.getValue();
				}
			else if(lcl.equals("interactorRef")) {
				p.interactorId = r.getElementText();
				}
			else if(lcl.equals("role")) {
				p.role = r.getElementText();
				}
			}
		else if(evt.isEndElement()) {
			final String lcl  = evt.asEndElement().getName().getLocalPart();
			if(lcl.equals("proteinParticipant")) break;
			}
		}
	if(p.interactorId==null) return null;
	return p;
	}



private static List<Participant> parseParticipantList(XMLEventReader r) throws XMLStreamException {
	final List<Participant> L = new ArrayList<>();
	while(r.hasNext()) {
		final XMLEvent evt = r.nextEvent();
		if(evt.isStartElement()) {
			final StartElement startE = evt.asStartElement();
			final String lcl  = startE.getName().getLocalPart();
			if(lcl.equals("proteinParticipant")) {
				final Participant p = parseParticipant(r);
				if(p!=null) L.add(p);
				}
			}
		else if(evt.isEndElement()) {
			final String lcl  = evt.asEndElement().getName().getLocalPart();
			if(lcl.equals("participantList")) break;
			}
		}
	return L;
	}

public static Interaction parse(XMLEventReader r, final StartElement root) throws XMLStreamException {
	final Interaction interaction = new Interaction();
	Attribute att = root.getAttributeByName(AbstractPsiObject.QNAME_ID);
	interaction.id = att.getValue();
	
	while(r.hasNext()) {
		final XMLEvent evt = r.nextEvent();
		if(evt.isStartElement()) {
			final StartElement startE = evt.asStartElement();
			final String lcl  = startE.getName().getLocalPart();
			if(lcl.equals("names")) {
				interaction.names = parseNames(r);
				}
			else if(lcl.equals("xref")) {
				interaction.references.addAll(parseXref(r));
				}
			else if(lcl.equals("participantList")) {
				interaction.participants.addAll(parseParticipantList(r));
				}
			}
		else if(evt.isEndElement()) {
			final String lcl  = evt.asEndElement().getName().getLocalPart();
			if(lcl.equals("interaction")) break;
			}
		}
	return interaction;
	}
}
