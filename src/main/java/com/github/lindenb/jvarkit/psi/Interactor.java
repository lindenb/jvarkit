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



import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;



public class Interactor extends AbstractPsiObject {
	
	public static class Organism  {
		String ncbiTaxId;
		Names names;
		}
	
	private Organism organism = null;
	
	private Interactor() {
	}

	
	public Organism getOrganism() {
		return this.organism;
	}
	
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Interactor)) return false;
		return this.getId().equals(Interactor.class.cast(obj).getId());
		}
	
	
	private static Organism paresOrganism(XMLEventReader r, final StartElement root) throws XMLStreamException {
		final Organism organism = new Organism();
		Attribute att = root.getAttributeByName(new QName("ncbiTaxId"));
		if(att!=null) organism.ncbiTaxId = att.getValue();
		while(r.hasNext()) {
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement()) {
				final String lcl  = evt.asStartElement().getName().getLocalPart();
				if(lcl.equals("names")) {
					organism.names = parseNames(r);
					}
				}
			else if(evt.isEndElement()) {
				final String lcl  = evt.asEndElement().getName().getLocalPart();
				if(lcl.equals("organism")) break;
				}
			}
		return organism;
		}
	
	
	protected static Interactor parseInteractor(final XMLEventReader r, final StartElement root) throws XMLStreamException {
		final Interactor interactor = new Interactor();
		Attribute att = root.getAttributeByName(AbstractPsiObject.QNAME_ID);
		interactor.id = att.getValue();
		
		while(r.hasNext()) {
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement startE = evt.asStartElement();
				final String lcl  = startE.getName().getLocalPart();
				if(lcl.equals("names")) {
					interactor.names = parseNames(r);
					}
				else if(lcl.equals("xref")) {
					interactor.references.addAll(parseXref(r));
					}
				else if(lcl.equals("organism")) {
					interactor.organism = paresOrganism(r, startE);
					}
				}
			else if(evt.isEndElement()) {
				final String lcl  = evt.asEndElement().getName().getLocalPart();
				if(lcl.equals("proteinInteractor")) break;
				}
			}
		return interactor;
		}

	
	}
