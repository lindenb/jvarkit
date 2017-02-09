/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

import htsjdk.variant.variantcontext.VariantContext;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.chart.Chart;
import javafx.scene.chart.PieChart;

public class VariantTypeChartFactory  extends VariantContextChartFactory
	{
	
		private static interface Category
			{
			public String getName();
			public boolean accept(final VariantContext ctx);
			}
		private  final Category categories[]=new Category[]{
			new Category(){
				@Override public String getName() {return "SNP Unfiltered";};
				@Override  public boolean accept(final VariantContext ctx) { return ctx.isSNP() && !ctx.isFiltered();}
				},
			new Category(){
				@Override public String getName() {return "SNP Filtered";};
				@Override  public boolean accept(final VariantContext ctx) { return ctx.isSNP() && ctx.isFiltered();}
				},
			new Category(){
				@Override public String getName() {return "Indel Unfiltered";};
				@Override  public boolean accept(final VariantContext ctx) { return ctx.isIndel() && !ctx.isFiltered();}
				},
			new Category(){
				@Override public String getName() {return "Indel Filtered";};
				@Override  public boolean accept(final VariantContext ctx) { return ctx.isIndel() && ctx.isFiltered();}
				},
			/* always LAST */
			new Category(){
				@Override public String getName() {return "Others";};
				@Override  public boolean accept(final VariantContext ctx) { return true;}
				}
			};
		private final long counts[]=new long[categories.length];
		
		@Override
		public String getName() { return "Variant Types";}
		
		@Override
		public void visit(VariantContext ctx)
			{
			for(int i=0;i< categories.length;++i)
				{
				if(categories[i].accept(ctx)) {
					counts[i]++;
					break;
					}
				}
			}
		
		@Override
		public Chart build() {	
	        final ObservableList<PieChart.Data> pieChartData = FXCollections.observableArrayList();
	        for(int i=0;i< categories.length;++i)
				{
	        	pieChartData.add( new PieChart.Data(categories[i].getName() +" "+counts[i], counts[i]));
				}
	        final PieChart chart = new PieChart(pieChartData);
	        chart.setTitle(this.getName());
	        return chart;
			}
	}
