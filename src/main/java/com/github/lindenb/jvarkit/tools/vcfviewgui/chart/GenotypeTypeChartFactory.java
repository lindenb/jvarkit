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

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class GenotypeTypeChartFactory
	extends VariantContextChartFactory {
	private final Map<String,Counter<GenotypeType>> sample2count=new TreeMap<>();
	@Override
	public String getName() {
		return "Genotype Type";
		}
	@Override
	public void visit(final VariantContext ctx) {
		for(final Genotype g:ctx.getGenotypes())
			{
			Counter<GenotypeType> count = this.sample2count.get(g.getSampleName());
			if(count==null) {
				count = new Counter<>();
				this.sample2count.put(g.getSampleName(), count);
				}
			count.incr(g.getType());
			}
		}
	@Override
	public Chart build() {	        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Sample");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        final List<XYChart.Series<String, Number>> gtype2count=new ArrayList<>(GenotypeType.values().length);
        for(final GenotypeType genotypeType :GenotypeType.values()) {
        	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
        	serie.setName(genotypeType.name());
        	gtype2count.add(serie);
        	
	        for(final String sampleName : this.sample2count.keySet())
	        	{
	        	serie.getData().add(new XYChart.Data<String,Number>(
	        			sampleName,
	        			this.sample2count.get(sampleName).count(genotypeType))
	        			);
	        	}
	        }
        
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle(this.getName());
        sbc.getData().addAll(gtype2count);
        sbc.setCategoryGap(0.2);
        return sbc;
        }
	
}
