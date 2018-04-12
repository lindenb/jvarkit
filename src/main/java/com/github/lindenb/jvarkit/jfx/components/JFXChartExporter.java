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
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;
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

	
private <T> String vectorR(
	final Iterable<T> data
	) {
	final StringBuilder sb = new StringBuilder("c(");
	boolean first=true;
	for(final T V : data){
		if(!first) sb.append(',');
		sb.append(V);
		first=false;
		};
	sb.append(")");
	return sb.toString();
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
	pw.print(vectorR(data.stream().map(D->D.getPieValue()).collect(Collectors.toList())));
	pw.print(",col=rainbow("+data.size()+")");
	pw.print(",labels=");
	pw.print(vectorR(data.stream().map(D->D.getName()).map(quoteR).collect(Collectors.toList())));
	title(chart);
	pw.println(")");
	}
private void exportBarChartSNToR(final BarChart<String, Number> chart) {
	final List<XYChart.Series<String,Number>> series = chart.getData();
	
	pw.print("barplot(table(");
	for(XYChart.Series<String,Number> serie : series) {
		pw.print(vectorR(serie.getData().stream().map(D->D.getYValue()).collect(Collectors.toList())));
		pw.print("),");
		}
	lab('x',chart.getXAxis());
	lab('y',chart.getYAxis());
	title(chart);
	pw.println(")");
	pw.flush();
	}

private void exportBarChartNSToR(final BarChart<Number, String> chart) {
}

private void exportStackedBarChartSNToR(final StackedBarChart<String, Number> chart) {
}

private void exportStackedBarChartNSToR(final StackedBarChart<Number, String> chart) {
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
	else if(chart instanceof StackedBarChart) {
		final StackedBarChart<?, ?> bc = StackedBarChart.class.cast(chart);
		if(bc.getData().isEmpty()) {
			LOG.warn("data is empty");
			return ;
			}
		if( bc.getXAxis() instanceof CategoryAxis &&
			bc.getYAxis() instanceof NumberAxis)
			{
			exportStackedBarChartSNToR((StackedBarChart<String, Number>)bc);
			return;
			}
		else if( bc.getXAxis() instanceof NumberAxis &&
			bc.getYAxis() instanceof CategoryAxis)
			{
			exportStackedBarChartNSToR((StackedBarChart<Number, String>)bc);
			return;
			}
		}
	else if(chart instanceof BarChart) {
		final BarChart<?, ?> bc = BarChart.class.cast(chart);
		if(bc.getData().isEmpty()) {
			LOG.warn("data is empty");
			return ;
			}
		if( bc.getXAxis() instanceof CategoryAxis &&
			bc.getYAxis() instanceof NumberAxis)
			{
			exportBarChartSNToR((BarChart<String, Number>)bc);
			return;
			}
		else if( bc.getXAxis() instanceof NumberAxis &&
			bc.getYAxis() instanceof CategoryAxis)
			{
			exportBarChartNSToR((BarChart<Number, String>)bc);
			return;
			}
		}
	LOG.warn("I don't know how to save a "+chart.getClass()+". Export is still experimental.");
	}
}
