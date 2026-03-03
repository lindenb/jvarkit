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
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.MiniList;
import com.google.gson.stream.JsonWriter;

/**
 * Scatter XY
 * plot one or multiple series of points(X,Y)
 * 
 */
public class ScatterXY extends  Chart implements MiniList<SeriesXY>  {
	private final List<SeriesXY> seriesXY;
	private final Double[] x_limits = new Double[] {null,null};
	private final Double[] y_limits = new Double[] {null,null};

	public ScatterXY() {
		this.seriesXY = new ArrayList<>();
		}
	
	public ScatterXY(List<SeriesXY> seriesXY) {
		this.seriesXY = new ArrayList<>(seriesXY);
		}
	public ScatterXY( SeriesXY seriesXY) {
		this(Collections.singletonList(seriesXY));
		}
	
	
	public ScatterXY setXLimits(final Double m,Double M) {
		this.x_limits[0] = m;
		this.x_limits[1] = M;
		return this;
		}
	public ScatterXY setYLimits(final Double m,Double M) {
		this.y_limits[0] = m;
		this.y_limits[1] = M;
		return this;
		}
	@Override
	public final SeriesXY get(int index) {
		return this.seriesXY.get(index);
		}
	@Override
	public final int size() {
		return this.seriesXY.size();
		}
	
	@Override
	public void saveXML(XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeStartElement("scatter");
		for(SeriesXY xy: this) {
			xy.saveXml(w);
			}
		w.writeEndElement();
		}
	
	
	private void savePlotlyHtml(XMLStreamWriter w) throws IOException,XMLStreamException {
		final String id= "plotid";
		w.writeStartDocument("UTF-8", "1.0");
		w.writeStartElement("html");
		w.writeStartElement("head");
		w.writeStartElement("script");
		//w.writeAttribute("src",this.plotly_library_url);
		w.writeEndElement();//script
		w.writeStartElement("script");
		w.writeAttribute("src","");
		w.writeCharacters("");
		w.writeEndElement();//script
		w.writeEndElement();//head
		w.writeStartElement("body");
		w.writeStartElement("div");
		w.writeAttribute("id", id);
		w.writeEndElement();///div
		w.writeEndElement();//body
		w.writeEndElement();//html
		w.writeEndDocument();
		w.close();
		}

	
	private void savePlotlyJS(Appendable w) throws IOException {
		w.append("var data").append(getId()).append(" = [")
			.append(this.stream().map(SERIES->SERIES.getId()).collect(Collectors.joining(",")))
			.append("];");
		w.append("var layout").append(getId()).append(" = {")
			.append(" title: {text: ")
			.append(StringUtils.doubleQuote(getTitle()))
			.append("}};");
		
		w.append("Plotly.newPlot('div")
			.append(getId())
			.append("', data")
			.append(getId())
			.append(", layout")
			.append(getId())
			.append(");");
		}
	
	public void savePlotly(Path directory,final String baseName) throws IOException {
		}
	
	public void saveMultiQC(final JsonWriter w) throws IOException {
		w.beginObject();
		
		w.name("data"); w.value(getId());
		w.name("plot_type"); w.value("linegraph");
		
		w.name("pconfig");
		w.beginObject();
			w.name("id"); w.value(getId());
			w.name("title"); w.value(getTitle());
			w.name("xlab"); w.value("xlab");
			w.name("ylab"); w.value("ylab");
		w.endObject();
		
		w.name("data");
		w.beginObject();
		for(SeriesXY series:this) {
			series.saveMultiQC(w);
			}
		w.endObject();//data
		w.endObject();
		}
	
	public void saveMultiQC(final Path filename) throws IOException {
		if(filename.getFileName().endsWith("_mqc.json")) {
			throw new IllegalArgumentException("filename should end with _mqc.json but got "+filename);
			}
		try(Writer w = Files.newBufferedWriter(filename)) {
			JsonWriter jw = new JsonWriter(w);
			saveMultiQC(jw);
			w.flush();
			}
		}
	
	
	
	}


