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
package com.github.lindenb.jvarkit.tools.simpleplot;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
@Program(
	name="simpleplot",
	description="Simple plots",
	keywords={"plot","chart"},
	creationDate = "20260220",
	modificationDate  = "20260220",
	jvarkit_hidden = true
	)
public class SimplePlot extends Launcher {
	private enum OutputMode {R,plotly,xml,tt};
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--list"},description="list factories and exit",help = true)
	private boolean list_factories =false;
	@Parameter(names={"-F"},description="Factory name")
	private String factoryName="undefined";
	@Parameter(names={"-m"},description="output mode")
	private OutputMode outputMode = OutputMode.R;
	@Parameter(names={"--plotly-url"},description="plotly javascript url")
	private String plotly_library_url = "https://cdn.plot.ly/plotly-3.3.0.min.js";
	
	@DynamicParameter(names = "-D", description = "Other parameters.")
	public Map<String, String> dynamicParams = new HashMap<String,String>() {{{
		put("xlab","__XLAB__");
		put("ylab","__YLAB__");
		put("main","__MAIN__");
		put("sub","__SUB__");
		}}};
	

	private static final Logger LOG = Logger.of(SimplePlot.class);
	private interface ChartFactory {
		public String getName();
		public default String getDescription() { return getName();}
		public int apply(BufferedReader br,PrintWriter out) throws IOException,XMLStreamException;
		}
	
	private abstract class AbstractChartFactory implements ChartFactory {
		private Function<String, List<String>> tokenizer = S->{
			return Collections.emptyList();
			};
		protected List<String> splitLine(final String s) {
			return tokenizer.apply(s);
			}
		protected Map.Entry<String,Double> splitSortUniq(final String s) {
			int i0=0;
			while(i0< s.length()) {
				if(!Character.isWhitespace(s.charAt(i0))) {
					break;
					}
				i0++;
				}
			int i1 = i0;
			while(i1< s.length()) {
				if(!Character.isDigit(s.charAt(i1))) {
					break;
					}
				i1++;
				}
			return new java.util.AbstractMap.SimpleEntry<String,Double>(null,null);
			}
		protected XMLStreamWriter openXMLStreamWriter(PrintWriter w) {
			return null;//TODO
			}
		
		public String getName() {
			return this.getClass().getSimpleName().replaceAll("Factory", "");
			}
		}
	
	private abstract class AbstractBoxPlot extends AbstractChartFactory {
		private class DataPoint {
			double value;
			String name;
			}
		public int apply(BufferedReader br,PrintWriter out)  throws IOException{
			String line;
			long n=1;
			final Map<String,List<DataPoint>> cat2points= new HashMap<>();
			while((line=br.readLine())!=null) {
				splitLine(line);
				String cat = "";
				List<DataPoint> L = cat2points.get(cat);
				if(L==null) {
					L=new ArrayList<>();
					cat2points.put(cat, L);
					}
				final DataPoint pt = new DataPoint();
				pt.value = 0.0;
				pt.name = "$"+n;
				L.add(pt);
				n++;
				}
			out.print("boxplot(c(");
			cat2points.entrySet().stream().map(P->String.valueOf(P.getValue())).collect(Collectors.joining(","));
			out.print("),names.arg=c(");
			cat2points.entrySet().stream().map(P->StringUtils.doubleQuote(P.getKey())).collect(Collectors.joining(","));
			out.print("))");
			return 0;
			}
		}

	
	private class SortUniqFactory extends AbstractChartFactory {
		
		public int apply(BufferedReader br,PrintWriter out)  throws IOException,XMLStreamException {
			String line;
			final List<Map.Entry<String,Double> > pairs = new ArrayList<>();
			while((line=br.readLine())!=null) {
				final Map.Entry<String,Double> pair = splitSortUniq(line);
				pairs.add(pair);
				}
			double max_v = pairs.stream().mapToDouble(P->P.getValue()).max().orElse(1.0);
			
			
			if(pairs.isEmpty()) {
				LOG.error("no data");
				return -1;
				}
			else if(SimplePlot.this.outputMode.equals(OutputMode.R)) {
				out.print("barplot(c(");
				pairs.stream().map(P->String.valueOf(P.getValue())).collect(Collectors.joining(","));
				out.print("),names.arg=c(");
				pairs.stream().map(P->StringUtils.doubleQuote(P.getKey())).collect(Collectors.joining(","));
				out.println("))");
				}
			else if(SimplePlot.this.outputMode.equals(OutputMode.plotly)) {
				final String id= "plotid";
				XMLStreamWriter w=openXMLStreamWriter(out);
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("html");
				w.writeStartElement("head");
				w.writeStartElement("script");
				w.writeAttribute("src", SimplePlot.this.plotly_library_url);
				w.writeEndElement();//script
				w.writeStartElement("script");
				w.writeCharacters(
						"var data"+id+" = [{type=\"bar\",x=["
								+ pairs.stream().map(P->StringUtils.doubleQuote(P.getKey())).collect(Collectors.joining(","))
								+ "],y=["
								+ pairs.stream().map(P->String.valueOf(P.getValue())).collect(Collectors.joining(","))
								+ "}]; Plotly.newPlot(" + StringUtils.doubleQuote(id)+", data"+id+"); "
						);
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
			else if(SimplePlot.this.outputMode.equals(OutputMode.tt)) {
				int name_len = pairs.stream().mapToInt(KV->KV.getKey().length()).max().orElse(0);
				
				for(Map.Entry<String,Double> p:pairs) {
					out.print(StringUtils.repeat(name_len - p.getKey().length(), ' '));
					out.print(p.getKey());
					out.print(" | ");
					out.print(String.format("%9ld",p.getValue().longValue()));
					out.print(" | ");
					out.println();
					}
				}
			else if(SimplePlot.this.outputMode.equals(OutputMode.xml)) {
				XMLStreamWriter w=openXMLStreamWriter(out);
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("barplot");
				for(Map.Entry<String,Double> p:pairs) {
					w.writeEmptyElement("bar");
					w.writeAttribute("name", p.getKey());
					w.writeAttribute("value", String.valueOf(p.getValue()));
					}
				w.writeEndElement();
				w.writeEndElement();
				w.close();
				}
			else
				{
				throw new IllegalArgumentException("not supported "+SimplePlot.this.outputMode);
				}	
			return 0;
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final List<ChartFactory> factories = new ArrayList<SimplePlot.ChartFactory>();
			if(list_factories) {
				for(ChartFactory f :factories) {
					stdout().println(f.getName()+"  |  "+f.getDescription());
					}
				return 0;
				}
			final ChartFactory factory  =factories.stream().filter(F->F.getName().equalsIgnoreCase(factoryName)).findFirst().orElse(null);
			if(factory==null) {
				LOG.error("Cannot find chart factory "+ factoryName);
				return -1;
				}
			int ret;
			try(BufferedReader br= super.openBufferedReader(super.oneFileOrNull(args))) {
				try(PrintWriter w= super.openPathOrStdoutAsPrintWriter(outputFile)) {
					ret = factory.apply(br, w);
					w.flush();
					}
				}
			return ret;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args) {
		new SimplePlot().instanceMainWithExit(args);
		}
	}
