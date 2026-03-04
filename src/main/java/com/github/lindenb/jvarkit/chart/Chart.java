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
import java.io.StringWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


public class Chart {
	private static int ID_GENERATOR=0;
	private final String id = String.valueOf("chart"+(++ID_GENERATOR));
	private String plotly_library_url = "https://cdn.plot.ly/plotly-3.3.0.min.js";

	private String mainTitle="";
	private String subTitle="";

	protected Chart() {
		}
	
	public String getId() {
		return id;
		}
	public String getTitle() {
		return mainTitle;
		}
	public Chart setTitle(final String title) {
		this.mainTitle = title;
		return this;
		}
	public Chart setSubTitle(final String title) {
		this.subTitle = title;
		return this;
		}
	
	public void setPlotlyLibraryUrl(String plotly_library_url) {
		this.plotly_library_url = plotly_library_url;
		}
	public String getPlotlyLibraryUrl() {
		return plotly_library_url;
		}
	public void saveXml(XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeComment("not implemented");
		}
	public void saveXml(Path p) throws IOException,XMLStreamException {
		final Charset charset=Charset.defaultCharset();
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		try(Writer w= Files.newBufferedWriter(p, charset)) {
			final XMLStreamWriter xw = xof.createXMLStreamWriter(w);
			xw.writeStartDocument(charset.displayName(), "1.0");
			saveXml(xw);
			xw.writeEndDocument();
			xw.flush();
			w.flush();
			}
		}
	
	public void savePlotlyJS(Appendable w) throws IOException {
		w.append("/* not implemented */");
		}
	
	public void savePlotly(final XMLStreamWriter w,final Charset charset) throws IOException,XMLStreamException {
		w.writeStartElement("html");
		w.writeStartElement("head");
		w.writeStartElement("title");
		w.writeCharacters(getTitle());
		w.writeEndElement();
		
		w.writeEmptyElement("meta");
		w.writeAttribute("charset", charset.displayName());
		w.writeEmptyElement("meta");
		w.writeAttribute("author", "Pierre Lindenbaum");
		w.writeEmptyElement("meta");
		w.writeAttribute("date", "todo");
	
		w.writeStartElement("script");
		w.writeAttribute("src",this.getPlotlyLibraryUrl());
		w.writeCharacters("");//empty
		w.writeEndElement();//script
		
		w.writeEndElement();//head
		w.writeStartElement("body");
		w.writeStartElement("h2");
		w.writeCharacters(getTitle());
		w.writeEndElement();//h2
		w.writeStartElement("div");
		w.writeAttribute("id", "div"+getId());
		w.writeEndElement();///div
		w.writeEmptyElement("hr");
		w.writeEndElement();//body
		
		w.writeStartElement("foot");
		
		w.writeStartElement("script");
		
		try(StringWriter sw=new StringWriter()) {
			savePlotlyJS(sw);
			w.writeCharacters(sw.toString());
			}
		w.writeEndElement();//script
		
		w.writeEndElement();//foot
		
		w.writeEndElement();//html
		}
	
	public void savePlotly(Path filename) throws IOException,XMLStreamException {
		if(filename.getFileName().endsWith(".html")) {
			throw new IllegalArgumentException("filename should end with .html but got "+filename);
			}
		final Charset charset=Charset.defaultCharset();
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		try(Writer w= Files.newBufferedWriter(filename, charset)) {
			final XMLStreamWriter xw = xof.createXMLStreamWriter(w);
			savePlotly(xw,charset);
			xw.flush();
			w.flush();
			}
		}
	}
