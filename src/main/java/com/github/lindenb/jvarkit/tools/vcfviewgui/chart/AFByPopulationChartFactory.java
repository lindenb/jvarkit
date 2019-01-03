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

import com.github.lindenb.jvarkit.tools.burden.MafCalculator;
import com.github.lindenb.jvarkit.tools.vcfviewgui.PedFile;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class AFByPopulationChartFactory extends VariantContextChartFactory {
	private static class Total
		{
		long num_maf=0L;
		double sum=0;
		}
	private final List<Total> popcount=new ArrayList<>(PedFile.Status.values().length);

	public AFByPopulationChartFactory()
		{
		for(int i=0;i< PedFile.Status.values().length;++i)
			{
			popcount.add(new Total());
			}
		}
	
	@Override
	public String getName() {
		return "Allele Frequency ( Population )";
		}
	@Override
	public void visit(final VariantContext ctx) {
		
		for(final Allele alt:ctx.getAlternateAlleles())
			{
			final MafCalculator mafCalculators[]=new MafCalculator[PedFile.Status.values().length];
			for(int i=0;i< mafCalculators.length;++i)
				{
				mafCalculators[i]=new MafCalculator(alt, ctx.getContig());
				}
			for(final Genotype gt:ctx.getGenotypes())
				{
				final PedFile.Sample sample=getPedigree().get(gt.getSampleName());
				if(sample==null)
					{
					mafCalculators[PedFile.Status.Unknown.ordinal()].
						add(gt,false);
					}
				else
					{
					mafCalculators[sample.getStatus().ordinal()].
						add(gt,sample.isMale());
					}
				}
			for(int i=0;i< mafCalculators.length;++i)
				{
				if(mafCalculators[i].isEmpty()) continue;
				final Total total=popcount.get(i);
				total.num_maf++;
				total.sum+=mafCalculators[i].getMaf();
				}
			}
		
		}

	@Override
	public Chart build() {	        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Population");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Mean-MAF");

    	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
    	serie.setName("Population");

        for(int i=0;i< this.popcount.size();++i) {
        	final PedFile.Status status = PedFile.Status.values()[i];
        	double v=0;
        	if(this.popcount.get(i).num_maf>0)
        		{
        		v=this.popcount.get(i).sum/((double)this.popcount.get(i).num_maf);
        		}
        	serie.getData().add(new XYChart.Data<String,Number>(
        			status.name(),v)
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
