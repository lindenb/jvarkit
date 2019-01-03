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

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class SamFlagsChartFactory extends ReadChartFactory
	{
    /** count position / base / number */
    private final Counter<Integer> flags2count= new Counter<>();
    
    @Override
    public String getName() {
    	return "Sam Flags";
    	}
    
    @Override
    public void visit(final SAMRecord rec) {
         _visit(rec.getFlags());	
        }
    @Override
    public void visit(final FastqRecord rec) {
         _visit(SAMFlag.READ_UNMAPPED.intValue());	
        }

   private void _visit(int flg) {
       	this.flags2count.incr(flg);
    	}
    
   @Override
    public StackedBarChart<String, Number> build()
        {	        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Flags");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        
    	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
    	serie.setName("Flags");
    	
    	
        for(final Integer L  : new TreeSet<Integer>(this.flags2count.keySet()))
        	{
        	serie.getData().add(new XYChart.Data<String,Number>(
        			String.valueOf(L),
        			this.flags2count.count(L))
        			);
        	}
        
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle("SAM Flags");
        sbc.getData().add(serie);
        sbc.setCategoryGap(0.2);
        sbc.setLegendVisible(false);
        return sbc;
        }
    }
