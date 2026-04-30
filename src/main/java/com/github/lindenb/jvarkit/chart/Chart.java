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
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.google.gson.Gson;
import com.google.gson.JsonElement;

/**
 * Base class for BarChart, ScatterXYChart , etc...
 */
public abstract class Chart {
	protected static int ID_GENERATOR=0;
	private String _id = null;
	private String plotly_library_url = "https://cdn.plot.ly/plotly-3.3.0.min.js";

	private String mainTitle="";
	private String subTitle="";

	protected Chart() {
		}
	
	public String getId() {
		if(StringUtils.isBlank(this._id)) {
			final String suffix = getTitle();
			this._id = String.valueOf("chart"+(++ID_GENERATOR))+(StringUtils.isBlank(suffix)?"":"."+suffix.replaceAll("[^A-Za-z0-9]+","_"));
			}
		return _id;
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
	public String getSubTitle() {
		return subTitle;
		}
	
	public void setPlotlyLibraryUrl(final String plotly_library_url) {
		this.plotly_library_url = plotly_library_url;
		}
	public String getPlotlyLibraryUrl() {
		return plotly_library_url;
		}
	public void saveXml(final XMLStreamWriter w) throws IOException,XMLStreamException {
		w.writeComment("not implemented");
		}
	public void saveXml(final Path p) throws IOException,XMLStreamException {
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
		if(!filename.getFileName().toString().endsWith(".html")) {
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
	public void saveR(final Appendable w) throws IOException {
		w.append("# not implemented\n");
		}
	
	public void saveR(final Path filename) throws IOException {
		try(PrintWriter pw = IOUtils.openPathForPrintWriter(filename)) {
			saveR(pw);
			pw.flush();
			}
		}
	
	protected abstract JsonElement getMultiQCJSon();
	
	public void saveMultiQC(final Path filename) throws IOException {
		if(!filename.getFileName().toString().endsWith("_mqc.json")) {
			throw new IllegalArgumentException("filename should end with _mqc.json but got "+filename);
			}
		final Gson gson = new Gson();
		final JsonElement root = getMultiQCJSon();
		try(PrintWriter pw= IOUtils.openPathForPrintWriter(filename)) {
			gson.toJson(root, pw);
			pw.flush();
			}
		}

	}
