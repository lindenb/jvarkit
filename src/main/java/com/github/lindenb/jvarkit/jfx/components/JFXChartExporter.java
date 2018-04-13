/*
The MIT License (MIT)
Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.jfx.components;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import com.github.lindenb.jvarkit.util.log.Logger;

import javafx.scene.chart.Axis;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.PieChart;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

/**
 * Experimentation class exporting JFX chart to R
 *
 */
public class JFXChartExporter {
private static final Logger LOG=Logger.build(JFXChartExporter.class).make();

	
private static Function<String, String> quoteR = (S)->{
	if(S==null) return "\"N/A\"";
	final StringBuilder sb = new StringBuilder(S.length()+2);
	sb.append("\"");
	for(int i=0;i< S.length();i++)
		{
		final char c = S.charAt(i);
		switch(c)
			{
			case '\'': sb.append("\\'");break;
			case '\"': sb.append("\"");break;
			case '\n': sb.append("\\n");break;
			case '\t': sb.append("\\t");break;
			case '\\': sb.append("\\\\");break;
			default: sb.append(c);break;
			}
		}
	sb.append("\"");
	return sb.toString();
	};
private final PrintWriter pw;


public JFXChartExporter(final PrintWriter pw)
	{
	this.pw = pw;
	}

private <T> void vectorR(
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
private void exportPieChartToR(final PieChart chart) {
	final List<PieChart.Data> data = chart.getData();
	pw.print("pie(");
	vectorR(data.stream().map(D->D.getPieValue()).iterator());
	pw.print(",col=rainbow("+data.size()+")");
	pw.print(",labels=");
	vectorR(data.stream().map(D->D.getName()).map(quoteR).iterator());
	title(chart);
	pw.println(")");
	}
private void exportBarChartSNToR(
		final XYChart<String, Number> chart,
		final boolean stacked
		) {
	final List<XYChart.Series<String,Number>> series = chart.getData();
	if(series.isEmpty()) {
		LOG.warn("empty chart");
		return;
	}
	pw.println("par(las="+(chart.getXAxis().getTickLabelRotation()==0?"0":"2")+")");
	pw.print("barplot(matrix(");
	
	vectorR(series.stream().
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
		vectorR(cataxis.getCategories().stream().map(quoteR).iterator());
		}
	if(series.size()>1) {
		pw.print(",legend=");
		vectorR(series.stream().map(T->T.getName()).map(quoteR).iterator());
		}
	pw.print(",beside="+(stacked?"F":"T"));
	lab('x',chart.getXAxis());
	lab('y',chart.getYAxis());
	title(chart);
	pw.println(")");
	pw.flush();
	}

private void exportBarChartNSToR(final XYChart<Number, String> chart,boolean stacked)  {

}


private void lab(char yx,final Axis<?> axis) {
	pw.print(",");
	pw.print(yx);
	pw.print("lab=");
	pw.print(quoteR.apply(axis.getLabel()));
	}

private void title(final Chart c) {
	pw.print(",main=");
	pw.print(quoteR.apply(c.getTitle()));
	}

@SuppressWarnings("unchecked")
/** export given chart to a R script, if possible */
public void exportToR(final Chart chart) {
	if(chart==null) return;
	if(chart instanceof PieChart) {
		exportPieChartToR(PieChart.class.cast(chart));
		return;
		}
	else if(chart instanceof StackedBarChart || chart instanceof BarChart) {
		final XYChart<?, ?> xyc = XYChart.class.cast(chart);
		final boolean stacked = (xyc instanceof StackedBarChart);
		if(xyc.getData().isEmpty()) {
			LOG.warn("data is empty");
			return ;
			}
		if( xyc.getXAxis() instanceof CategoryAxis &&
			xyc.getYAxis() instanceof NumberAxis)
			{
			exportBarChartSNToR((XYChart<String, Number>)xyc,stacked);
			return;
			}
		else if( xyc.getXAxis() instanceof NumberAxis &&
				xyc.getYAxis() instanceof CategoryAxis)
			{
			exportBarChartNSToR((XYChart<Number, String>)xyc,stacked);
			return;
			}
		}
	LOG.warn("I don't know how to save a "+chart.getClass()+". Export is still experimental.");
	}
}
