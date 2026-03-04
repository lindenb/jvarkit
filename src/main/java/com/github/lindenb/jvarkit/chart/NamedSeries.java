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
package com.github.lindenb.jvarkit.chart;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

/** a bar in a bar plot */
public class NamedSeries extends AbstractSeries<NamedY> {
	public NamedSeries(final String name, List<NamedY> values) {
		super(name,values);
		}
	public NamedSeries(final String name,double value) {
		super(name,Collections.singletonList(new NamedY("data",value)));
		}
	public double sum() {
		return stream().mapToDouble(NamedY::getY).sum();
		}
	public NamedY getNamedYByName(final String s) {
		return stream().filter(NY->NY.getName().equals(s)).findFirst().orElse(null);
		}
	void saveXml(XMLStreamWriter w) throws XMLStreamException {
		w.writeStartElement("series");
		w.writeAttribute("name", this.getName());
		w.writeAttribute("id", this.getId());
		w.writeAttribute("size",String.valueOf(this.size()));
		for(NamedY ny: this) {
			ny.saveXml(w);
			}
		w.writeEndElement();
		}
	@Override
	public String toString() {
		return getName() + "=["+stream().map(X->X.toString()).collect(Collectors.joining(","))+"]";
		}
	}
