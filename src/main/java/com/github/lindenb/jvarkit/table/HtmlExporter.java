/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.table;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.net.UrlSupplier;


public class HtmlExporter extends AbstractTableExporter {
	private UrlSupplier urlSupplier = new UrlSupplier();
	
	public UrlSupplier getUrlSupplier() {
		return urlSupplier;
		}
	
	public void setUrlSupplier(final UrlSupplier urlSupplier) {
		this.urlSupplier = urlSupplier;
		}
	
	
	public void print(final Table t,final XMLStreamWriter w) throws XMLStreamException
		{
		w.writeStartElement("html");
		w.writeStartElement("body");
		w.writeStartElement("table");
		w.writeStartElement("thead");
		w.writeStartElement("caption");
		w.writeCharacters(""+t.getTitle());
		w.writeEndElement();
		w.writeStartElement("tr");
		for(int x=0;x< t.getColumnCount();++x)
			{
			final String str = t.getColumn(x).getName();
			w.writeStartElement("th");
			w.writeCharacters(str);
			w.writeEndElement();
			}
		w.writeEndElement();//end tr
		w.writeEndElement();//thead
		w.writeStartElement("tbody");
		for(int y=0;y< t.getRowCount();++y)
			{
			final List<Object> row= t.getRow(y);
			w.writeStartElement("tr");
			for(int x=0;x< row.size();++x)
				{
				w.writeStartElement("td");
				final String str = t.getColumn(x).toString(row.get(x));
				final UrlSupplier.LabelledUrl url=getUrlSupplier()==null?null:getUrlSupplier().of(str).stream().findFirst().orElse(null);
				if(url!=null)
					{
					w.writeStartElement("a");
					w.writeAttribute("title", url.getLabel());
					w.writeAttribute("href", url.getUrl());
					}
				w.writeCharacters(str);
				if(url!=null)
					{
					w.writeEndElement();
					}
				
				w.writeEndElement();//tr
				}
			w.writeEndElement();//tr
			}
		w.writeEndElement();//tbody
		w.writeEndElement();//table
		
		w.writeEndElement();//body
		w.writeEndElement();//html
		}
	
	@Override
	public  void saveTableTo(final Table table,final PrintWriter pw) throws IOException {
		final XMLOutputFactory xof = XMLOutputFactory.newInstance();
		try {
			final XMLStreamWriter w=xof.createXMLStreamWriter(pw);
			this.print(table,w);
			w.flush();
			w.close();
			}
		catch(final XMLStreamException err) {
			throw new IOException(err);
			}
		}
}
