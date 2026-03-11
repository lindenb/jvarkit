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
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.MiniList;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

public class BarPlot extends AbstractChartXY  implements MiniList<NamedSeries>  {
	private final List<NamedSeries> delegate;
	private boolean beside=false;
	
	public BarPlot(final NamedSeries series) {
		this(Collections.singletonList(series));
		}
	public BarPlot(List<NamedSeries> L) {
		this.delegate = new ArrayList<>(L);
		}
	public BarPlot(final Counter<String> counter) {
		this(counter.stream()
				.map(KV->new NamedSeries(KV.getKey(), KV.getValue().doubleValue()))
				.collect(Collectors.toList())
			);
		}
	
	public void sortOnName() {
		final SmartComparator cmp = new SmartComparator();
		Collections.sort(delegate,(A,B)->cmp.compare(A.getName(), B.getName()));
		}
	
	@Override
	public NamedSeries get(int index) {
		return this.delegate.get(index);
		}
	@Override
	public int size() {
		return  this.delegate.size();
		}
	public boolean isBeside() {
		return beside;
		}
	public void setBeside(boolean beside) {
		this.beside = beside;
		}
	
	private List<String> getCategories() {
		return new ArrayList<>( this.stream()
			.flatMap(NS->NS.stream())
			.map(NS->NS.getName())
			.collect(Collectors.toCollection(LinkedHashSet::new))
			);
		}
	
	
	@Override
	public void saveR(Appendable w) throws IOException {
		// TODO Auto-generated method stub
		super.saveR(w);
		}
	
	
	

	@Override
	public JsonObject 	getMultiQCJSon()
		{
		final  JsonObject o = new JsonObject();
		o.add("id", new JsonPrimitive(getId()));
		o.add("plot_type", new JsonPrimitive("bargraph"));
		
		final  JsonObject pconfig = new JsonObject();
		o.add("pconfig", pconfig);
		
		pconfig.add("id", new JsonPrimitive(getId()));
		pconfig.add("title", new JsonPrimitive(getTitle()));
		pconfig.add("xlab", new JsonPrimitive(getXAxisLabel()));
		pconfig.add("ylab", new JsonPrimitive(getYAxisLabel()));
		//pconfig.add("xlog", new JsonPrimitive(isLogX()));
		pconfig.add("ylog", new JsonPrimitive(isLogY()));
		
		final List<String> categories = getCategories();
		final  JsonObject data = new JsonObject();
		o.add("data", data);
		for(NamedSeries series:this) {
			JsonObject datao = new JsonObject();
			data.add(series.getName(), datao);
			for(final String cat: categories) {
				final NamedY y = series.getNamedYByName(cat);
				if(y==null) {
					datao.addProperty(cat, 0.0);
					}
				else
					{
					datao.addProperty(cat, y.getY());
					}
				}
			}
		return o;
		}
	
	
	@Override
	public void saveXml(final XMLStreamWriter w) throws IOException, XMLStreamException {
		w.writeStartElement("barplot");
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
		w.writeStartElement("categories");
		w.writeAttribute("count", String.valueOf(getCategories().size()));
		for(String cat: getCategories()) {
			w.writeStartElement("category");
			w.writeCharacters(cat);
			w.writeEndElement();
			}
		w.writeEndElement();
		w.writeStartElement("series");
		w.writeAttribute("count", String.valueOf(this.size()));
		for(NamedSeries ns: this) {
			ns.saveXml(w);
			}
		w.writeEndElement();
		w.writeEndElement();
		}
	
	public void savePlotlyJS(final Appendable w) throws IOException {
		List<String> categories_names = getCategories();
		if(categories_names.isEmpty()) return;
		
	
			for(int cat_index = 0; cat_index < categories_names.size();++cat_index) {
				w.append("trace").append(getId()+"_"+cat_index).append(" = {");
				w.append("x:[");
				for(int x=0;x<size();x++) {
					if(x>0) w.append(",");
					final NamedSeries ns = get(x);
					w.append(StringUtils.doubleQuote(ns.getName()));
					}
				w.append("], y:[");
				for(int x=0;x<size();x++) {
					if(x>0) w.append(",");
					final NamedSeries ns = get(x);
					final NamedY ny = ns.getNamedYByName(categories_names.get(cat_index));
					if(ny!=null) {
						w.append(String.valueOf(ny.getY()));
						}
					else
						{
						w.append("0");
						}
					}
				w.append("], name:");
				w.append(StringUtils.doubleQuote(categories_names.get(cat_index)));
				w.append(", type:'bar'");
				w.append("};\n");
				}
			
			w.append("var data").append(getId()).append("= [");
			for(int index=0;index< categories_names.size();++index) {
				if(index>0) w.append(",");
				w.append("trace"+getId()+"_"+index);
				}
			w.append("];\n");
			
			w.append("var layout"+getId()+" = {");
			if(!isBeside()) {
				w.append("barmode:"+StringUtils.doubleQuote("stack"));
				}
			w.append("};\n");
			
			w.append("Plotly.newPlot('div")
				.append(getId())
				.append("', data")
				.append(getId())
				.append(", layout")
				.append(getId())
				.append(");\n");
			
		
		}
	}
