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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class TiTvChartFactory  extends VariantContextChartFactory
	{
	private static class TiTv{
		int transition=0;
		int transvertion=0;
		private boolean isATGC(char c)
			{
			return c=='A' || c=='T' || c=='G' || c=='C';
			}
		void watch(final List<Allele> alleles)
			{
			if(alleles==null || alleles.size()!=2) return;
			char ref='\0';
			char alt='\0';
			for(Allele a: alleles)
				{
				if(a.isNoCall()) return;
				if(a.isSymbolic()) return;
				if(!a.isCalled()) return;
				String display=a.getDisplayString().toUpperCase();
				if(display.length()!=1) return ;
				char c=display.charAt(0);
				if(!isATGC(c)) return;
				if(a.isReference())
					{
					if(ref!='\0') return;
					ref=c;
					}
				else
					{
					if(alt!='\0') return;
					alt=c;
					}
				}
			if(ref=='\0' || alt=='\0') return;
			
			
	        if ((ref == 'A' && alt == 'G')
	                || (ref == 'G' && alt == 'A')
	                || (ref == 'C' && alt == 'T')
	                || (ref == 'T' && alt == 'C')
	                )
	        	{
	        	transition++;
	        	}
	        else
	        	{
	        	transvertion++;
	        	}
	    
			}
		};
	
	private final Map<String,TiTv> sample2titv=new TreeMap<>();
	private final TiTv all=new TiTv();
	
	@Override
	public String getName() {
		return "Ti/Tv";
		}
	@Override
	public void visit(final VariantContext ctx) {
		if(this.sample2titv.isEmpty() && ctx.getNSamples()>0)
			{
			for(final String s:ctx.getSampleNames())
				{
				this.sample2titv.put(s, new TiTv());
				}
			}
		for(final Genotype g:ctx.getGenotypes())
			{
			this.sample2titv.get(g.getSampleName()).watch(g.getAlleles());
			}
		this.all.watch(ctx.getAlleles());
		}
	@Override
	public Chart build() {	        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Sample");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        final List<XYChart.Series<String, Number>> type2count =new ArrayList<>(3);
        for(int i=0;i<2;++i) {
        	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
        	serie.setName(i==0?"Transition":"Transversion");
        	type2count.add(serie);
        	
        	if( this.sample2titv.isEmpty())
	        	{
	        	serie.getData().add(new XYChart.Data<String,Number>(
	        			"ALL",
	        			(i==0?all.transition:all.transvertion)
	        			));
	        	}
        	else
	        	{
		        for(final String sampleName : this.sample2titv.keySet())
		        	{
		        	final TiTv titv= this.sample2titv.get(sampleName);
		        	serie.getData().add(new XYChart.Data<String,Number>(
		        			sampleName,
		        			(i==0?titv.transition:titv.transvertion)
		        			));
		        	}
		        }
        	}
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle(this.getName());
        sbc.getData().addAll(type2count);
        sbc.setCategoryGap(0.2);
        return sbc;
        }
	
}
