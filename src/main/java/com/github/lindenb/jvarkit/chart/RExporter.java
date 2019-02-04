package com.github.lindenb.jvarkit.chart;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
 import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

public class RExporter extends ChartExporter {
	private static final Logger LOG=Logger.build(RExporter.class).make();

		
	private static Function<String, String> quoteR = (S)->{
		if(S==null) return "\"N/A\"";
		final StringBuilder sb = new StringBuilder(S.length()+2);
		sb.append("\"");
		sb.append(StringUtils.escapeC(S));
		sb.append("\"");
		return sb.toString();
		};

	private <T> void vectorR(
		final PrintWriter pw,
		final Iterator<T> it
		) {
		pw.print("c(");
		boolean first=true;
		while(it.hasNext()){
			if(!first) pw.append(',');
			pw.print(it.next());
			first=false;
			};
		pw.print(")");
		}

	private String colors(int n) {
		return  "rainbow("+n+")";
		}

	/**

	 barplot(
		 matrix(
		c(10,5,15,25,11,12,13,14),
		ncol=4,byrow=TRUE
		),
		names.arg=c("A","B","C","D"),
			legend = c("X1","X2"),
		main="T",xlab="xxx", ylab="yy", beside=TRUE)
	 
	 */
	private void exportPieChartToR(
			final PrintWriter pw,
			final PieChart chart
			) {
		final List<PieChart.Data> data = chart.getData();
		if(data.isEmpty()) {
			LOG.warn("empty chart");
			return;
			}
		pw.print("pie(");
		vectorR(pw,data.stream().map(D->D.getPieValue()).iterator());
		pw.print(",col="+colors(data.size()));
		pw.print(",labels=");
		vectorR(pw,data.stream().map(D->D.getName()).map(quoteR).iterator());
		title(pw,chart);
		pw.println(")");
		}

	private void exportScatterNNToR(
			final PrintWriter pw,
			final XYChart<Number, Number> chart) 
		{
		final List<XYChart.Series<Number,Number>> series = chart.getData().stream().
				filter(P->!P.getData().isEmpty()).
				collect(Collectors.toList());
		if(series.isEmpty()) {
			LOG.warn("empty chart");
			return;
			}
		par(pw,chart.getXAxis());
		int serie_index=0;
		final char type=(chart instanceof ScatterChart?'p':'o');
		while(serie_index < series.size())
			{
			if(serie_index>0) {
				pw.println("par(new=TRUE)");
				}
			
			pw.print("plot(matrix(");
			vectorR(pw,series.get(serie_index).getData().
					stream().
					map(D->String.valueOf(D.getXValue())+","+String.valueOf(D.getYValue())).
					iterator()
					);
			pw.print(",ncol=2,byrow=TRUE)");
			if(serie_index>0) {
				pw.print(",axes=FALSE");
				}
			lim(pw,'x',NumberAxis.class.cast(chart.getXAxis()),chart,D->D.getY());
			lim(pw,'y',NumberAxis.class.cast(chart.getYAxis()),chart,D->D.getY());
			pw.print(",pch="+serie_index);
			pw.print(",type=\""+type+"\"");
			lab(pw,'x',chart.getXAxis());
			lab(pw,'y',chart.getYAxis());
			if(serie_index==0)
				{
				title(pw,chart);
				}
			if(series.size()>1) {
				pw.print(",col="+colors(series.size()));
				}
			pw.println(")");
			
			serie_index++;
			}
		if(series.size()>1) {
			pw.print("legend(\"topright\",legend=");
			vectorR(pw,series.stream().map(T->T.getName()).map(quoteR).iterator());
			pw.print(",pch=");
			vectorR(pw,IntStream.range(0, series.size()).mapToObj(T->Integer.toString(T)).iterator());
			pw.println(")");
			}
			
		pw.flush();
		}

	private void exportBarChartSNToR(
			final PrintWriter pw,
			final XYChart<String, Number> chart
			
			) {
		final boolean stacked = chart instanceof StackedBarChart;
		final List<XYChart.Series<String,Number>> series = chart.getData();
		if(series.isEmpty()) {
			LOG.warn("empty chart");
			return;
		}
		par(pw,chart.getXAxis());
		pw.print("barplot(matrix(");
		
		vectorR(pw,series.stream().
				flatMap(S->S.getData().stream()).
				map(D->D.getYValue()).
				iterator()
				);
		pw.print(",ncol=");
		pw.print(series.get(0).getData().size());
		pw.print(",byrow=TRUE)");
		if(chart.getXAxis() instanceof CategoryAxis)
			{
			final CategoryAxis cataxis = CategoryAxis.class.cast(chart.getXAxis());
			pw.print(",names.arg=");
			vectorR(pw,cataxis.getCategories().stream().map(quoteR).iterator());
			}
		if(series.size()>1) {
			pw.print(",legend=");
			vectorR(pw,series.stream().map(T->T.getName()).map(quoteR).iterator());
			}
		if(series.size()>1) {
			pw.print(",col="+colors(series.size()));
			}
		
		pw.print(",beside="+(stacked?"F":"T"));
		lab(pw,'x',chart.getXAxis());
		lab(pw,'y',chart.getYAxis());
		title(pw,chart);
		pw.println(")");
		pw.flush();
		}


	private <T extends Number > void exportHeatMap(final PrintWriter pw,
		final HeatMapChart<T> map) {
		if(map.isEmpty()) {
			LOG.warn("empty chart");
			return;
			}
		final List<String> xlabels = map.getXAxis().getCategories();
		final List<String> ylabels =  map.getYAxis().getCategories();
		if(xlabels.size()<2 || ylabels.size()<2) {
			LOG.warn("HeatMap must have at least 2 rows and 2 columns");
			return;
			}
		boolean first=true;
		pw.print("heatmap(matrix(c(");
		for(final String y: ylabels) {
			for(final String x: xlabels) {
					final Optional<T> v=map.get(x, y);
					if(!first) pw.print(",");
					first=false;
					pw.print(v.isPresent()?String.valueOf(v.get()):"NA");
				}
			}
		pw.print("),ncol=" +xlabels.size()+",nrow="+ylabels.size());
		pw.print(",dimnames=list(");
		vectorR(pw, ylabels.stream().map(quoteR).iterator());
		pw.print(",");
		vectorR(pw, xlabels.stream().map(quoteR).iterator());
		pw.print(")),scale=\"none\",Colv=NA,Rowv=NA");
		lab(pw,'x', map.getXAxis());
		lab(pw,'y', map.getYAxis());
		title(pw, map);
		pw.println(")");
		}

	private void exportBarChartNSToR(final PrintWriter pw,final XYChart<Number, String> chart)  {
		LOG.warn("I don't know how to save a "+chart.getClass()+". Export is still experimental.");
	}


	private void lab(final PrintWriter pw,char yx,final Axis<?> axis) {
		pw.print(",");
		pw.print(yx);
		pw.print("lab=");
		pw.print(quoteR.apply(axis.getLabel()));
		}
	private void lim(final PrintWriter pw,char yx,final NumberAxis axis,XYChart<?, ?> chart,final Function<XYChart.Data<?, ?>,Object> extractor) {
		pw.print(",");
		pw.print(yx);
		pw.print("lim=c(");
		pw.print(lowerBound(axis,chart,extractor));
		pw.print(",");
		pw.print(upperBound(axis,chart,extractor));
		pw.print(")");
		}
	private void title(final PrintWriter pw,final Chart c) {
		pw.print(",main=");
		pw.print(quoteR.apply(c.getTitle()));
		}

	private void par(final PrintWriter pw,final Axis<?> xaxis) {
		pw.println("par(las="+(xaxis.getTickLabelRotation()==0?"0":"2")+")");
		}

	@SuppressWarnings("unchecked")
	/** export given chart to a R script, if possible */
	public void exportToR(final PrintWriter pw,final Chart chart) {
		if(chart==null) return;
		chart.update();
		if(chart instanceof PieChart) {
			exportPieChartToR(pw,PieChart.class.cast(chart));
			return;
			}
		else if(chart instanceof HeatMapChart) {
			final HeatMapChart<? extends Number> heat = HeatMapChart.class.cast(chart);
			exportHeatMap(pw, heat);
			return;
			}
		else if(chart instanceof ScatterChart ||
				chart instanceof LineChart)
			{
			final XYChart<?, ?> xyc = XYChart.class.cast(chart);
			if(xyc.getYAxis().isNumber() && xyc.getYAxis().isNumber())
				{
				exportScatterNNToR(pw,(XYChart<Number, Number>)xyc);
				return;
				}
			}
		else if(chart instanceof StackedBarChart || chart instanceof BarChart) {
			final XYChart<?, ?> xyc = XYChart.class.cast(chart);
			if(xyc.getData().isEmpty()) {
				LOG.warn("data is empty");
				return ;
				}
			if(xyc.getXAxis().isString() &&
				xyc.getYAxis().isNumber())
				{
				exportBarChartSNToR(pw,(XYChart<String, Number>)xyc);
				return;
				}
			else if( xyc.getXAxis().isNumber() &&
					xyc.getYAxis().isString())
				{
				exportBarChartNSToR(pw,(XYChart<Number, String>)xyc);
				return;
				}
			}
		LOG.warn("I don't know how to save a "+chart.getClass()+". Export is still experimental.");
		}

	private double lowerBound(final NumberAxis axis,final XYChart<?, ?> chart,final Function<XYChart.Data<?, ?>,Object> extractor) {
		if(axis.getLowerBound()!=null) return axis.getLowerBound();
		return chart.getData().stream().
			flatMap(S->S.getData().stream()).
			map(extractor).
			mapToDouble(O->Number.class.cast(O).doubleValue()).
			min().orElse(0.0);
		}
	private double upperBound(final NumberAxis axis,final XYChart<?, ?> chart,final Function<XYChart.Data<?, ?>,Object> extractor) {
		if(axis.getUpperBound()!=null) return axis.getUpperBound();
		return chart.getData().stream().
			flatMap(S->S.getData().stream()).
			map(extractor).
			mapToDouble(O->Number.class.cast(O).doubleValue()).
			max().orElse(0.0);
		}

	
}
