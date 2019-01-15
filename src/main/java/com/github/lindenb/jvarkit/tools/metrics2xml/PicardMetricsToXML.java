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

package com.github.lindenb.jvarkit.tools.metrics2xml;

import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
/**
## Example  

beware the **Locale**: For example I generated the file using a french Locale. I had to transform the comma to dot.

```bash
$ cat jeter.metrics |\
   tr "," "." |\
   java -jar dist/metrics2xml.jar |\
   xsltproc ./src/main/resources/xsl/picardmetrics2html.xsl - 

<html xmlns:p="http://picard.sourceforge.net/">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<script type="text/javascript" src="https://www.google.com/jsapi"></script><scri
 google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChartidp24896);
      function drawChartidp24896() {
        var dataidp24896 = google.visualization.arrayToDataTable([
          ['insert_size','All_Reads.fr_count','All_Reads.rf_count'],[2,3979.0,38
        ]);

        var optidp24896 = {
          title: 'stdin'
        };

        var chartidp24896 = new google.visualization.LineChart(document.getEleme
        chartidp24896.draw(dataidp24896, optidp24896);
        }
</script>
</head>
<body><div>
<h1>stdin</h1>
<h2>Headers</h2>
<dl>
<dt>net.sf.picard.metrics.StringHeader</dt>
<dd>net.sf.picard.analysis.CollectMultipleMetrics INPUT=../align/sorted.bam REFE
<dt>net.sf.picard.metrics.StringHeader</dt>
<dd>Started on: Tue Nov 19 11:31:58 CET 2013</dd>
</dl>
<h2>net.sf.picard.analysis.InsertSizeMetrics</h2>
<table border="1">
<thead>
<th>MEDIAN_INSERT_SIZE</th>
<th>MEDIAN_ABSOLUTE_DEVIATION</th>
<th>MIN_INSERT_SIZE</th>
<th>MAX_INSERT_SIZE</th>
<th>MEAN_INSERT_SIZE</th>
<th>STANDARD_DEVIATION</th>
<th>READ_PAIRS</th>
<th>PAIR_ORIENTATION</th>
<th>WIDTH_OF_10_PERCENT</th>
<th>WIDTH_OF_20_PERCENT</th>
<th>WIDTH_OF_30_PERCENT</th>
<th>WIDTH_OF_40_PERCENT</th>
<th>WIDTH_OF_50_PERCENT</th>
<th>WIDTH_OF_60_PERCENT</th>
<th>WIDTH_OF_70_PERCENT</th>
<th>WIDTH_OF_80_PERCENT</th>
<th>WIDTH_OF_90_PERCENT</th>
<th>WIDTH_OF_99_PERCENT</th>
(...)
```

### An XML-to-HTML example

![https://jvarkit.googlecode.com/svn/wiki/picardmetrics01.jpg](https://jvarkit.googlecode.com/svn/wiki/picardmetrics01.jpg)

 */
@Program(name="picardmetrics2xml",
	description="transforms a picard metrics file to XML. See http://plindenbaum.blogspot.fr/2013/02/making-use-of-picard-metrics-files.html",
	keywords={"picard","xml","metrics"}
	)
public class PicardMetricsToXML
	extends Launcher
	{
	private static final Logger LOG=Logger.build(PicardMetricsToXML.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output=null;
	@Parameter(names="-s",description="print sum")
	private boolean print_sums=false;


	private static final String NS="http://picard.sourceforge.net/";
	private static class MinMax
		{
		Double min=null;
		Double max=null;
		double sum=0.0;
		int count=0;
		}
	
	private MetricsFile<MetricBase, Comparable<?>> metricsFile;
	private XMLStreamWriter out;
	
	private MetricsFile<MetricBase, Comparable<?>> getFile() {
		return metricsFile;
		}
	
	private void printHeaders() throws XMLStreamException
		{
		out.writeStartElement("headers");
		for (final Header h : getFile().getHeaders())
			{
			out.writeStartElement("header");
            out.writeAttribute("class", h.getClass().getName());
            out.writeCharacters(h.toString());
            out.writeEndElement();
			}
		out.writeEndElement();
		}
	
	private void printBeanMetrics() throws XMLStreamException
		{
		out.writeStartElement("metrics");
		if(getFile().getMetrics()!=null && !getFile().getMetrics().isEmpty())
			{
			Class<?> beanType=getFile().getMetrics().get(0).getClass();
			Field fields[]=beanType.getFields();
			List<MinMax> minmaxList=new ArrayList<MinMax>(fields.length);
			for(int i=0;i< fields.length;++i) minmaxList.add(new MinMax());
			out.writeStartElement("thead");
			out.writeAttribute("class",beanType.getName());
			for(Field field:fields )
				{
				out.writeStartElement("th");
				out.writeAttribute("class",field.getType().getName());
				out.writeCharacters(field.getName());
				out.writeEndElement();//th
				}
			out.writeEndElement();//thead
			out.writeStartElement("tbody");
			for(MetricBase bean:getFile().getMetrics())
				{
				out.writeStartElement("tr");
				for(int i=0; i<fields.length; ++i)
					{
					try{
						final Object value = fields[i].get(bean);
						
						if(value==null)
							{
							out.writeEmptyElement("td");
							out.writeAttribute("xsi",XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI,"nil","true");
							}
						else
							{
							if(print_sums && value instanceof Number)
								{
								double numeric=Number.class.cast(value).doubleValue();
								MinMax minmax=minmaxList.get(i);
								if(!(Double.isNaN(numeric) || Double.isInfinite(numeric)))
									{
									if(minmax.count==0)
										{
										minmax.max=numeric;
										minmax.min=numeric;										}
									else
										{
										minmax.max=Math.max(numeric,minmax.max);
										minmax.min=Math.min(numeric,minmax.min);
										}
									minmax.sum+=numeric;
									minmax.count++;
									}
								}
							out.writeStartElement("td");
							out.writeCharacters(String.valueOf(value));
							out.writeEndElement();//td
							}
						}
					catch (IllegalAccessException e)
						{
						out.writeEmptyElement("td");
						out.writeAttribute("error",e.getClass().getSimpleName());
						}
					
					}
				out.writeEndElement();//tr
				}
			out.writeEndElement();//tbody
			
			if(print_sums)
				{
				out.writeStartElement("tfoot");
				for(int i=0; i<fields.length; ++i)
					{
					MinMax minmax=minmaxList.get(i);
					out.writeEmptyElement("measurement");
					out.writeAttribute("for", fields[i].getName());
					if(minmax.count>0)
						{
						out.writeAttribute("count", String.valueOf(minmax.count));
						out.writeAttribute("min", String.valueOf(minmax.min));
						out.writeAttribute("max", String.valueOf(minmax.max));
						out.writeAttribute("mean", String.valueOf(minmax.sum/minmax.count));
						}		
					}
				out.writeEndElement();//tfoot
				}
			}	
		out.writeEndElement();
		}
	@SuppressWarnings({"unchecked","rawtypes"})
	private void printHistogram()throws XMLStreamException
		{
		
		final List<Histogram<?>> nonEmptyHistograms = new ArrayList<Histogram<?>>();
		for(Histogram<?> histo: getFile().getAllHistograms())
			{
			if(histo!=null && !histo.isEmpty()) nonEmptyHistograms.add(histo);
			}
		if(nonEmptyHistograms.isEmpty()) return;
		SortedSet keys = new TreeSet(nonEmptyHistograms.get(0).comparator());
		for(Histogram<?> histo:nonEmptyHistograms)
			{
			keys.addAll(histo.keySet());
			}
		out.writeStartElement("histogram");
		out.writeAttribute("class",  nonEmptyHistograms.get(0).keySet().iterator().next().getClass().getName());
		out.writeStartElement("thead");
		out.writeStartElement("th");
		out.writeCharacters(String.valueOf(nonEmptyHistograms.get(0).getBinLabel()));
		out.writeEndElement();//th
		for (final Histogram<?> histo : nonEmptyHistograms)
		 	{
			out.writeStartElement("th");
			out.writeCharacters(String.valueOf(histo.getValueLabel()));
			out.writeEndElement();//th
		 	}
		out.writeEndElement();//thead
		out.writeStartElement("tbody");
		for(final Object key: keys)
			{
			out.writeStartElement("tr");
			out.writeStartElement("td");
			out.writeCharacters(String.valueOf(key));
			out.writeEndElement();//td
			for (final Histogram histo : nonEmptyHistograms)
				{
				Histogram.Bin<? extends Comparable> bin = (Histogram.Bin<? extends Comparable>)histo.get((Comparable)key);
				if(bin==null)
					{
					out.writeEmptyElement("td");
					out.writeAttribute("xsi",XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI,"nil","true");
					}
				else
					{
					out.writeStartElement("td");
					out.writeCharacters(String.valueOf(bin.getValue()));
					out.writeEndElement();//td
					}
				}
			out.writeEndElement();//tr
			}
		out.writeEndElement();
		out.writeEndElement();
		}
	
	private void parse(String filename,Reader in) throws IOException,XMLStreamException
		{
		LOG.info("processing "+filename);
		out.writeStartElement("metrics-file");
		out.writeAttribute("file", filename);
		this.metricsFile = new MetricsFile<MetricBase, Comparable<?>>();
		this.metricsFile.read(in);
		printHeaders();
		printBeanMetrics();
		printHistogram();
		out.writeEndElement();
		}
	
	private void parse(File file) throws IOException,XMLStreamException
		{
		parse(file.toString(),new FileReader(file));
		}
	@Override
	public int doWork(List<String> args)
		{
		OutputStream pstream=null;
		try
			{
			pstream = super.openFileOrStdoutAsStream(output);
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			this.out= xmlfactory.createXMLStreamWriter(pstream,"UTF-8");
			this.out.setDefaultNamespace(NS);
			out.writeStartDocument("UTF-8","1.0");
			out.writeStartElement("picard-metrics");
			out.writeDefaultNamespace(NS);
			out.writeAttribute(XMLConstants.XMLNS_ATTRIBUTE, XMLConstants.XML_NS_URI,"xsi", XMLConstants.W3C_XML_SCHEMA_INSTANCE_NS_URI);

			if(args.isEmpty())
				{
				parse("stdin",new InputStreamReader(stdin()));
				}
			else
				{
				for(final String filename:args)
					{
					parse(new File(filename));
					}
				}
			out.writeEndElement();
			out.writeEndDocument();
			out.flush();
			out.close();
			pstream.flush();
			pstream.close();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.out);
			CloserUtil.close(pstream);
			}
		}
	
	public static void main(String[] args)
		throws IOException,XMLStreamException
		{
		new PicardMetricsToXML().instanceMainWithExit(args);
		}
	}
