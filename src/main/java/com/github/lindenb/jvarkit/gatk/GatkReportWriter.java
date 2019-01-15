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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.gatk;

import java.io.PrintStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportDataType;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import htsjdk.samtools.util.RuntimeIOException;

public abstract class GatkReportWriter {
public static enum Format { DEFAULT,HTML,XML,TSV};

public static GatkReportWriter createWriter(final Format fmt) 
	{
	switch(fmt)
		{
		case HTML: return createHtmlWriter();
		case XML: return createXmlWriter();
		case TSV: return createTsvWriter();
		default: return createDefaultWriter();
		}
	}


public static GatkReportWriter createDefaultWriter()
	{
	return new GatkReportWriter() {
		@Override
		public void print(final GATKReport report,final PrintStream out) {
			report.print(out);
			out.flush();
			}
		};
	}

public static GatkReportWriter createHtmlWriter()
	{
	return new HtmlGatkReportWriter();
	}

public static GatkReportWriter createXmlWriter()
	{
	return new XmlGatkReportWriter();
	}


public abstract void print(final GATKReport report,final PrintStream out);


private static class HtmlGatkReportWriter extends GatkReportWriter
	{
	@Override
	public void print(final GATKReport report,final PrintStream out) {
		XMLStreamWriter w=null;
		try {
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w=xof.createXMLStreamWriter(out);
			w.writeStartElement("div");
			/* loop over each table */
			for(final GATKReportTable table:report.getTables()) {
				w.writeStartElement("div");
				w.writeStartElement("table");
				
				/* write table header */
				w.writeStartElement("thead");
				w.writeStartElement("caption");
				w.writeCharacters(table.getTableName());
				w.writeEmptyElement("br");
				w.writeCharacters(table.getTableDescription());
				w.writeEndElement();
				w.writeStartElement("tr");
				/* loop over each column */
				for(final GATKReportColumn col:table.getColumnInfo())
					{
					w.writeStartElement("th");
					w.writeCharacters(col.getColumnName());
					w.writeEndElement();//th
					}
				w.writeEndElement();//tr
				w.writeEndElement();//thread
				
				/* write table body */
				w.writeStartElement("tbody");
				/* loop over each row */
				for(int i=0;i< table.getNumRows();++i)
					{
					w.writeStartElement("tr");
					/* loop over each column */
					for(int j=0;j< table.getNumColumns();++j)
						{
						w.writeStartElement("td");
						w.writeCharacters(
								obj2str(
									table.getColumnInfo().get(j),
									table.get(i, j)
									)
								);
						w.writeEndElement();//th
						}
					w.writeEndElement();//tr
					}
				
				w.writeEndElement();//tbody
				w.writeEndElement();//table
				w.writeEndElement();//div
				}
			w.writeEndElement();//div
			w.flush();
			out.flush();
			}
		catch(final Exception err)
			{
			throw new RuntimeIOException(err);
			}
		}
	}

private static class XmlGatkReportWriter extends GatkReportWriter
{
@Override
public void print(final GATKReport report,final PrintStream out) {
	XMLStreamWriter w=null;
	try {
		XMLOutputFactory xof=XMLOutputFactory.newFactory();
		w=xof.createXMLStreamWriter(out);
		w.writeStartElement("report");
		/* loop over each table */
		for(final GATKReportTable table:report.getTables()) {
			w.writeStartElement("table");
			w.writeAttribute("name", table.getTableName());
			w.writeAttribute("description", table.getTableDescription());
			w.writeAttribute("rows",String.valueOf(table.getNumRows()));
			w.writeAttribute("columns",String.valueOf(table.getNumColumns()));
			
			
			w.writeStartElement("columns");
			/* loop over each column */
			for(int j=0;j< table.getNumColumns();++j)
				{
				final GATKReportColumn col = table.getColumnInfo().get(j);
				w.writeEmptyElement("column");
				w.writeAttribute("name", col.getColumnName());
				w.writeAttribute("index",String.valueOf(j));
				w.writeAttribute("format", String.valueOf(col.getFormat()));
				w.writeAttribute("type", String.valueOf(col.getDataType().name()));
				
				boolean nileable=false;
				for(int i=0;i< table.getNumRows();++i)
					{
					if(table.get(i, j)==null ) {nileable=true;break;}
					}
				w.writeAttribute("null", String.valueOf(nileable));
				
				Class<?> onlyClass=null;
				for(int i=0;i< table.getNumRows();++i)
					{
					final Object o = table.get(i, j);
					if(o==null )continue;
					if(onlyClass!=null && !onlyClass.equals(o.getClass()))
							{
							onlyClass=null;
							break;
							}
					onlyClass=o.getClass();
					}
				if(onlyClass!=null) w.writeAttribute("class", onlyClass.getName());
				
				boolean number=false;
				for(int i=0;i< table.getNumRows();++i)
					{
					final Object o = table.get(i, j);
					if(o==null )continue;
					if(!(o instanceof Number)) {number=false;break;}
					number=true;
					}
				w.writeAttribute("number", String.valueOf(number));
				}
			w.writeEndElement();//columns
			
			/* write table body */
			w.writeStartElement("rows");
			/* loop over each row */
			for(int i=0;i< table.getNumRows();++i)
				{
				w.writeStartElement("row");
				w.writeAttribute("index",String.valueOf(i));
				/* loop over each column */
				for(int j=0;j< table.getNumColumns();++j)
					{
					final Object obj = table.get(i, j);
					w.writeStartElement("cell");
					w.writeAttribute("column",String.valueOf(j));
					w.writeAttribute("class",obj==null?"null":obj.getClass().getName());					
					if(obj==null) w.writeAttribute("null","true");					
					w.writeCharacters(obj==null?"null":obj.toString());
					w.writeEndElement();//cell
					}
				w.writeEndElement();//row
				}
			
			w.writeEndElement();//rows
			w.writeEndElement();//table
			}
		w.writeEndElement();//report
		w.flush();
		out.flush();
		}
	catch(final Exception err)
		{
		throw new RuntimeIOException(err);
		}
	}
}

public static GatkReportWriter createTsvWriter()
{
return new GatkReportWriter()
		{
		@Override
		public void print(GATKReport report, PrintStream out) {

			boolean firstTable=true;
			/* loop over each table */
			for(final GATKReportTable table:report.getTables()) {
				if(!firstTable) out.println();
				firstTable=false;
				out.print("##Name\t");
				out.println(table.getTableName());
				out.print("##Description\t");
				out.println(table.getTableDescription());
				/* loop over each column */
				for(int j=0;j< table.getNumColumns();++j)
					{
					final GATKReportColumn col = table.getColumnInfo().get(j);
					out.print(j==0?"#":"\t");
					out.print(col.getColumnName());
					}
				out.println();
				/* loop over each row */
				for(int i=0;i< table.getNumRows();++i)
					{
					/* loop over each column */
					for(int j=0;j< table.getNumColumns();++j)
						{
						final Object obj = table.get(i, j);
						if(j>0) out.print("\t");
						if(obj!=null) out.print(obj.toString());
						}
					out.println();
					}
				
				}
			out.flush();			
			}
		};
}



private static String obj2str(final GATKReportColumn col,final Object obj) {
	 final String value;
	if ( obj == null )
        value =  "null";
    else if ( col.getDataType().equals(GATKReportDataType.Unknown) && (obj instanceof Double || obj instanceof Float) )
    	value =   String.format("%.8f", obj);
    else
    	value =    String.format(col.getFormat(), obj);
	
	return String.format(col.getColumnFormat().getValueFormat(), value);
	}
}
