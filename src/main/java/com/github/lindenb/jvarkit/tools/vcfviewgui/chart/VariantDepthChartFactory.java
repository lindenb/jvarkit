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
import htsjdk.variant.vcf.VCFConstants;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class VariantDepthChartFactory extends VariantContextChartFactory {
	private final Counter<Integer> afindexcount=new Counter<>();
	
	
	@Override
	public String getName() {
		return "Variant Depth";
		}
	@Override
	public void visit(final VariantContext ctx)
			{
			if(ctx.hasAttribute(VCFConstants.DEPTH_KEY))
				{
				final int v;
				try
				{
				v=ctx.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
				if(v<0) return;
				this.afindexcount.incr(v-v%10);
				}
			catch(NumberFormatException err) {
				return;
				}
			}
		else
			{
			//LOG.warning("NOT AF in variant");
			}
		}

	@Override
	public Chart build() {	        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("DP");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        
    	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
    	serie.setName(xAxis.getLabel());
    	
    	for(Integer dp: new TreeSet<>(this.afindexcount.keySet()))
        	{
        	serie.getData().add(new XYChart.Data<String,Number>(
        			String.valueOf(dp),
        			this.afindexcount.count(dp))
        			);
        	}
    	
        final BarChart<String, Number> sbc =
                new BarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle(this.getName());
        sbc.getData().add(serie);
        sbc.setCategoryGap(0.2);
        sbc.setLegendVisible(false);
        return sbc;
        }
	
}
