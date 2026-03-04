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
public class ScatterXY extends  AbstractChartXY implements MiniList<SeriesXY>  {
	private final List<SeriesXY> seriesXY;
	public ScatterXY() {
		this.seriesXY = new ArrayList<>();
		}
	
	public ScatterXY(List<SeriesXY> seriesXY) {
		this.seriesXY = new ArrayList<>(seriesXY);
		}
	public ScatterXY( SeriesXY seriesXY) {
		this(Collections.singletonList(seriesXY));
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
	public void saveXml(final XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeStartElement("scatter");
		w.writeAttribute("id", getId());
		
		w.writeStartElement("title");
		w.writeCharacters(getTitle());
		w.writeEndElement();
		w.writeStartElement("x-axis");
		w.writeAttribute("log", String.valueOf(isLogX()));
		w.writeCharacters(getXAxisLabel());
		w.writeEndElement();
		w.writeStartElement("y-axis");
		w.writeAttribute("log", String.valueOf(isLogY()));
		w.writeCharacters(getYAxisLabel());
		w.writeEndElement();		
		
		
		for(SeriesXY xy: this) {
			xy.saveXml(w);
			}
		w.writeEndElement();
		}
	
	/*
	 {
		    type: 'line',
		    xref: 'x',
		    yref: 'paper',   // "paper" makes it span the full plot height
		    x0: 2e6,
		    x1: 2e6,
		    y0: 0,
		    y1: 1,
		    line: { color: 'red', width: 1, dash: 'dot' }
		  }, */
	
	protected void writePlotlyShapes(Appendable w) throws IOException { 
		w.append("[]");
		}
	
	@Override
	public void savePlotlyJS(Appendable w) throws IOException {
		w.append("\n");
		for(SeriesXY series:this) {
			series.plottlyJS(w);
			}
		
		w.append("var data").append(getId()).append(" = [")
			.append(this.stream().map(SERIES->SERIES.getId()).collect(Collectors.joining(",")))
			.append("];\n");
		w.append("var layout").append(getId()).append(" = {")
			.append(" title: {text: ")
			.append(StringUtils.doubleQuote(getTitle()))
			.append("},")
			.append("xaxis: {")
			.append("title: {text: ").append(StringUtils.doubleQuote(getXAxisLabel())).append("}")
			.append(",autorange:true");
			if(super.isLogX()) {
				w.append(",type:'log'");
				}
		
		w.append("},")
			.append("yaxis: {")
			.append("title: {text: ").append(StringUtils.doubleQuote(getYAxisLabel())).append("}")
			.append(",autorange:true");
		if(super.isLogY()) {
			w.append(",type:'log'");
			}
		w.append("},shapes:");
		
		writePlotlyShapes(w);
		
		w.append("};\n");
		
		w.append("Plotly.newPlot('div")
			.append(getId())
			.append("', data")
			.append(getId())
			.append(", layout")
			.append(getId())
			.append(");\n");
		}
	
	
	
	public void saveMultiQC(final JsonWriter w) throws IOException {
		w.beginObject();
		
		w.name("data"); w.value(getId());
		w.name("plot_type"); w.value("linegraph");
		
		w.name("pconfig");
		w.beginObject();
			w.name("id"); w.value(getId());
			w.name("title"); w.value(getTitle());
			w.name("xlab"); w.value(getXAxisLabel());
			w.name("ylab"); w.value(getYAxisLabel());
		w.endObject();
		
		w.name("data");
		w.beginObject();
		for(SeriesXY series:this) {
			series.saveMultiQC(w);
			}
		w.endObject();//data
		w.endObject();
		}
	@Override
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


