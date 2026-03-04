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
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.MiniList;

public class BarPlot extends AbstractChartXY  implements MiniList<NamedSeries>  {
	private final List<NamedSeries> delegate;
	private List<String> categories = null;
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
	
	public void setCategoriesName(List<String> categories)  {
		this.categories = new ArrayList<String>(categories);
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
	
	private String getCategoryAt(int index) {
		if(this.categories==null || index<0 || index>=this.categories.size()) {
			return "#cat"+(1+index);
			}
		return this.categories.get(index);
		}
	
	public void savePlotlyJS(final Appendable w) throws IOException {
		final int n_categories = stream().mapToInt(L->L.size()).max().orElse(0);
		if(n_categories==0) return;
		
		if(n_categories>1) {
			for(int cat_index = 0; cat_index < n_categories;++cat_index) {
				
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
					if(cat_index< ns.size()) {
						w.append(String.valueOf(ns.get(cat_index)));
						}
					else
						{
						w.append("0");
						}
					}
				w.append("], name:");
				w.append(StringUtils.doubleQuote(getCategoryAt(cat_index)));
				w.append(", type:'bar'");
				w.append("};\n");
				}
			
			w.append("var data").append(getId()).append("= [");
			for(int index=0;index< n_categories ;++index) {
				if(index>0) w.append(",");
				w.append(StringUtils.doubleQuote("trace"+getId()+"_"+index));
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
		else
			{
			w.append("data").append(getId()).append(" = [{");
			w.append("x:[");
			w.append(this.stream().map(S->StringUtils.doubleQuote(S.getName())).collect(Collectors.joining(",")));
			w.append("], y:[");
			w.append(this.stream().map(S->String.valueOf(S.get(0))).collect(Collectors.joining(",")));
			w.append("], type:'bar'");
			w.append("}];\n");
			
			w.append("Plotly.newPlot('div")
				.append(getId())
				.append("', data")
				.append(getId())
				//.append(", layout")
				//.append(getId())
				.append(");\n");
			}
	}
	}
