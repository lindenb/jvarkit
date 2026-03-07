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

import java.awt.geom.Point2D;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.google.gson.JsonArray;
import com.google.gson.stream.JsonWriter;

public class DataXY extends AbstractDataY {
	double x;
	public DataXY(final double x, final double y) {
		super(y);
		this.x = x;
		}
	public DataXY(Point2D pt) {
		this(pt.getX(),pt.getY());
		}
	
	
	public double getX() {
		return x;
		}
	
	JsonArray buildMultiQCJson() {
		JsonArray a=new JsonArray();
		a.add(getX());
		a.add(getY());
		return a;
		}
	
	void saveMultiQC(final JsonWriter w) throws IOException {
		w.beginArray();
		w.value(getX());
		w.value(getY());
		w.endArray();
		}
	
	void saveXml(XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeEmptyElement("point");
		w.writeAttribute("x", String.valueOf(getX()));
		w.writeAttribute("y", String.valueOf(getY()));
		}
	}