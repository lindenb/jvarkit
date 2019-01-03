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
package com.github.lindenb.jvarkit.tools.vcfviewgui.chart;

import java.util.TreeSet;

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.variant.variantcontext.VariantContext;
import javafx.scene.chart.Chart;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class VariantQualChartFactory  extends VariantContextChartFactory
		{
		private final Counter<Integer> qual2count=new Counter<>();
		
		@Override
		public String getName() { return "Variant QUAL";}
		
		@Override
		public void visit(final VariantContext ctx)
			{
			if(ctx.hasLog10PError())
				{
				qual2count.incr((int)ctx.getPhredScaledQual());
				}
			else
				{
				qual2count.incr(-1);
				}
			}
		
		@Override
		public Chart build() {
	    	final NumberAxis xAxis = new NumberAxis();
	    	xAxis.setLabel("QUAL");
	        final NumberAxis yAxis = new NumberAxis();
	        yAxis.setLabel("Count");

	    	final  XYChart.Series<Number,Number> serie = new XYChart.Series<Number,Number>();
	    	serie.setName(xAxis.getLabel());
	    	
	        for(final Integer i: new TreeSet<>(this.qual2count.keySet()))
	        	{
	        	serie.getData().add(new XYChart.Data<Number,Number>(
	        			i,this.qual2count.count(i))
	        			);
	        	}
	        
	        
	        final LineChart<Number, Number> sbc =
	                new LineChart<Number, Number>(xAxis, yAxis);
	        sbc.setTitle(this.getName());
	        sbc.getData().add(serie);
	        sbc.setCreateSymbols(false);
	        sbc.setLegendVisible(false);
	        return sbc;
	        }
	}
