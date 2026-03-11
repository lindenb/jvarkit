package com.github.lindenb.jvarkit.chart;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.MiniList;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

public class BoxPlotChart extends AbstractChartXY implements MiniList<NamedSeries> {
	private final List<NamedSeries> delegate;
	public BoxPlotChart(final NamedSeries series) {
		this(Collections.singletonList(series));
		}
	public BoxPlotChart(List<NamedSeries> L) {
		this.delegate = new ArrayList<>(L);
		}
	public BoxPlotChart(final Counter<String> counter) {
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

	

	@Override
	public JsonObject 	getMultiQCJSon()
		{
		final  JsonObject o = new JsonObject();
		o.add("id", new JsonPrimitive(getId()));
		o.add("plot_type", new JsonPrimitive("boxplot"));
		
		final  JsonObject pconfig = new JsonObject();
		o.add("pconfig", pconfig);
		
		pconfig.add("id", new JsonPrimitive(getId()));
		pconfig.add("title", new JsonPrimitive(getTitle()));
		pconfig.add("xlab", new JsonPrimitive(getXAxisLabel()));
		pconfig.add("ylab", new JsonPrimitive(getYAxisLabel()));
		pconfig.add("xlog", new JsonPrimitive(isLogX()));
		pconfig.add("ylog", new JsonPrimitive(isLogY()));
		
		final  JsonObject data = new JsonObject();
		o.add("data", data);
		for(NamedSeries series:this) {
			data.add(series.getName(), series.getMultiQCArrayY());
			}
		return o;
		}
	
	
	
	
	public void savePlotlyJS(final Appendable w) throws IOException {
		if(isEmpty()) return;
		for(NamedSeries ns : this) {
			w.append("trace").append(ns.getId()).append(" = {");
			w.append("x:[");
			w.append(ns.stream().map(NY->String.valueOf(NY.getY())).collect(Collectors.joining(",")));
			w.append("], type:'box',  boxpoints: 'Outliers', name:");
			w.append(StringUtils.doubleQuote(ns.getName()));
			w.append("};\n");
			}
		
			
		w.append("var data").append(getId()).append("= [");
		w.append(stream().map(ns->"trace"+ns.getId()).collect(Collectors.joining(",")));
		w.append("];\n");
			
		w.append("var layout"+getId()+" = {");
		w.append("title:{ text: ").append(StringUtils.doubleQuote(getTitle())).append("}");
		
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
