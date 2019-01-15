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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class AlleleFrequencyChartFactory extends VariantContextChartFactory {
	private final Counter<Integer> afindexcount=new Counter<>();
	private final List<Limit> limits =new ArrayList<>();
	private int nSamples=0;
	private static class Limit
		{
		final int index;
		final double v;
		final String label;
		Limit(int index,double v,String label) {this.index=index;this.v=v;this.label=label;}
		}
	
	@Override
	public String getName() {
		return "Allele Frequency";
		}
	@Override
	public void visit(final VariantContext ctx) {
		if(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
			{
			if(this.limits.isEmpty())
    			{
				final DecimalFormat df = new DecimalFormat("#.###");
				int ndiv;
    			if(ctx.getNSamples()>0)
    				{
    				this.nSamples=ctx.getNSamples();
    				ndiv=(2*this.nSamples);
    				}
    			else
    				{
    				ndiv = 10;
    				}
    			
    			for(int x=0;x< ndiv;++x)
    				{
    				double v= x*(1.0/ndiv);
    				String cat="??";
    				if(x==0)
						{
						cat="<"+df.format(v);
						}
    				else if(x+1==ndiv)
						{
						cat=">="+df.format(v);
						}
					else
						{
						cat="["+df.format(v)+" - "+df.format(v+ (1.0/ndiv))+"[";
						}
					
    				this.limits.add(new Limit(this.limits.size(),v,cat));
    				}
    			}
			for(final Double v: super.getAttributeAsDoubleList(ctx,VCFConstants.ALLELE_FREQUENCY_KEY))
				{
				Limit cat=null;
				for(int x=0;x<limits.size();++x)
					{
					if(x==0 && v< limits.get(0).v)
						{
						cat = limits.get(0);
						break;
						}
					else if(x+1==limits.size())
						{
						cat= limits.get(x);
						break;
						}
					else if( v>=limits.get(x).v && v<limits.get(x+1).v)
						{
						cat=  limits.get(x);
						break;
						}
					
					}		
				if(cat!=null)
					{
					this.afindexcount.incr(cat.index);
					}
				else
					{
					//LOG.warning("?? AF: "+v);
					}
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
    	xAxis.setLabel("AF");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        
    	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
    	serie.setName(xAxis.getLabel());
    	
    	
        for(final Limit limit  :this.limits)
        	{
        	if(this.afindexcount.count(limit.index)==0L) continue;
        	serie.getData().add(new XYChart.Data<String,Number>(
        			limit.label,
        			this.afindexcount.count(limit.index))
        			);
        	}
        
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle("Allele Frequency"+(this.nSamples>0?" (Nbr. Sample(s):"+this.nSamples+")":""));
        sbc.getData().add(serie);
        sbc.setCategoryGap(0.2);
        sbc.setLegendVisible(false);
        return sbc;
        }
	
}
