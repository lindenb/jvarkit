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
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class BasesPerPositionChartFactory extends ReadChartFactory
	{
	/** all bases seen so far */
    private final Set<Character> all_chars=new TreeSet<>();
    /** count position / base / number */
    private final List<Counter<Character>> pos2base2count= new ArrayList<>();
    
    @Override
    public String getName() {
    	return "Bases/Position";
    	}
    
    @Override
    public void visit(final SAMRecord rec) {
    	byte bases[]= rec.getReadBases();
    	if(bases==null || bases.length==0) return ;
    	
    	if(! rec.getReadUnmappedFlag() && rec.getReadNegativeStrandFlag())
    		{
    		bases = Arrays.copyOf(bases, bases.length);//because it would modify the read itself.
    		SequenceUtil.reverseComplement(bases);
    		}
        _visit(bases);	
        }
    
    @Override
    public void visit(final FastqRecord rec) {
    	_visit(rec.getReadString().getBytes());
    	}
    
   private void _visit(final byte bases[]) {
        if(bases==null || bases.length==0) return;
        
        while(this.pos2base2count.size()<bases.length)
        	{
        	this.pos2base2count.add(new Counter<>());
        	}
        
    	for(int x=0;x< bases.length;++x)
    		{
    		char c=(char)bases[x];
    		this.pos2base2count.get(x).incr(c);
    		this.all_chars.add(c);
    		}
    	return;
    	}
    @Override
    public StackedBarChart<String, Number> build()
        {        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Position in Read");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        
        final List<XYChart.Series<String, Number>> base2count=new ArrayList<>(all_chars.size());
        for(final Character base:all_chars) {
        	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
        	serie.setName(base.toString());
        	base2count.add(serie);
        	
	        for(int i=0;i<  this.pos2base2count.size();++i)
	        	{
	        	serie.getData().add(new XYChart.Data<String,Number>(
	        			String.valueOf(i+1),
	        			this.pos2base2count.get(i).count(base))
	        			);
	        	}
	        }
        
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle("Position/Base/Count");
        sbc.getData().addAll(base2count);
        sbc.setCategoryGap(0.2);
        
        return sbc;
        }
    }
