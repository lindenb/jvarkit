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
		
		private static class CategoryImpl implements Category
		{
			final boolean filtered;
			final boolean snp;
			final boolean id;
			CategoryImpl(boolean filtered,boolean snp,boolean id) {
				this.filtered=filtered;
				this.snp=snp;
				this.id=id;
			}
		@Override public String getName()
			{
			return (snp?"SNP":"INDEL")+" "+
					(filtered?"":"NOT")+" FILTERed "+
					(id?"HAS":"NO")+" ID";
			}
		@Override public boolean accept(final VariantContext ctx) {
			if(snp && !ctx.isSNP()) return false;
			if(!snp && !ctx.isIndel()) return false;
			if(filtered!=ctx.isFiltered()) return false;
			if(id!=ctx.hasID()) return false;
			return true;
			}
		}

		
		private  final Category categories[]=new Category[]{
			new CategoryImpl(true,true,true),
			new CategoryImpl(true,true,false),
			new CategoryImpl(true,false,true),
			new CategoryImpl(true,false,false),
			new CategoryImpl(false,true,true),
			new CategoryImpl(false,true,false),
			new CategoryImpl(false,false,true),
			new CategoryImpl(false,false,false),
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
	        	if(counts[i]==0) continue;
	        	pieChartData.add( new PieChart.Data(categories[i].getName() +" "+counts[i], counts[i]));
				}
	        final PieChart chart = new PieChart(pieChartData);
	        chart.setTitle(this.getName());
	        chart.setLegendVisible(false);
	        return chart;
			}
	}
