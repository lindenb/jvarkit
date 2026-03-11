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

import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.google.gson.JsonArray;
import com.google.gson.stream.JsonWriter;

public  class SeriesXY extends AbstractSeries<DataXY> {
	//private Color color;
	public SeriesXY() {
		this("",new ArrayList<DataXY>());
		}
	
	public SeriesXY(String name) {
		this(name,Collections.emptyList());
		}
	
	public SeriesXY( List<DataXY> data) {
		this("",data);
		}
	public SeriesXY(String name, List<DataXY> data) {
		super(name,data);
		}
	
	
	
	OptionalDouble getMinX() {
		return stream().mapToDouble(DataXY::getX).min();
		}
	OptionalDouble getMaxX() {
		return stream().mapToDouble(DataXY::getX).max();
		}
	OptionalDouble getMinY() {
		return stream().mapToDouble(DataXY::getY).min();
		}
	OptionalDouble getMaxY() {
		return stream().mapToDouble(DataXY::getY).max();
		}
	
	
	
	void saveXml(XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeStartElement("series");
		if(!StringUtils.isBlank(getName())) {
			w.writeAttribute("name", getName());
			}
		w.writeAttribute("size", String.valueOf(size()));
		if(!isEmpty()) {
			w.writeAttribute("min-x", String.valueOf(getMinX().getAsDouble()));
			w.writeAttribute("max-x", String.valueOf(getMaxX().getAsDouble()));
			w.writeAttribute("min-y", String.valueOf(getMinY().getAsDouble()));
			w.writeAttribute("max-y", String.valueOf(getMaxY().getAsDouble()));
			for(DataXY p:this) {
				p.saveXml(w);
				}
			}
		
		
		
		w.writeEndElement();
		}
	void plottlyJS(Appendable w) throws IOException {
		w.append("var "+getId())
			.append(" = {\nx:[")
			.append(this.stream().map(PT->String.valueOf(PT.getX())).collect(Collectors.joining(",")))
			.append("],\ny:[")
			.append(this.stream().map(PT->String.valueOf(PT.getY())).collect(Collectors.joining(",")))
			.append("]");
		if(!StringUtils.isBlank(getName())) {
			w.append(",\nname:")
				.append(StringUtils.doubleQuote(getName()));
			}
			
		w.append("};\n");
		}
	
	/** reorder DATA on X and then Y */
	public SeriesXY sort() {
		Collections.sort(super.delegate,(A,B)->{
			final int i = Double.compare(A.getX(), B.getX());
			if(i!=0) return i;
			return Double.compare(A.getY(), B.getY());
			});
		return this;
		}
	public SeriesXY normalize() {
		final double minY = getMinY().orElse(0);
		final double maxY = getMaxY().orElse(0);
		double distance = maxY-minY;
		if(distance==0) distance=1;
		for(DataXY pt:this) {
			pt.setY((pt.getY()-minY)/distance);
			}
		return this;
		}
	
	Map.Entry<String,JsonArray> buildMultiQCJson() {
		final JsonArray array = new JsonArray();
		for(DataXY xy:this) {
			array.add(xy.buildMultiQCJson());
			}
		return new AbstractMap.SimpleEntry<>(this.getName(),array);
		}
	
	
	public void saveMultiQC(JsonWriter w) throws IOException {
		w.name(getName());
		w.beginArray();
		for(DataXY xy:this) {
			xy.saveMultiQC(w);
			}
		w.endArray();
		}
	}