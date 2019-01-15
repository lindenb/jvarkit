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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class QualityPerPositionChartFactory 
	extends ReadChartFactory{
	
	private  static class QualData {
		double sum=0.0;
		long count=0L;
		byte m=Byte.MAX_VALUE;
		byte M=Byte.MIN_VALUE;
		void visit(byte qual) {
			this.sum+=qual;
			this.count++;
			if(qual<m) m=qual;
			if(qual>M) M=qual;
			}
		double get() { return this.sum/this.count;}
		double min() { return m;}
		double max() { return M;}
		}
    /** count position / base / number */
    private final List<QualData> pos2qual= new ArrayList<>();
    
    
    @Override
    public String getName() {
    	return "Quality/Read Pos.";
    	}
    
    @Override
    public void visit(final SAMRecord rec) {
    	byte quals[]= rec.getBaseQualities();
    	if(quals==null || quals.length==0) return ;
    	
    	if(! rec.getReadUnmappedFlag() && rec.getReadNegativeStrandFlag())
    		{
    		quals = Arrays.copyOf(quals, quals.length);//because it would modify the read itself.
    		SequenceUtil.reverseQualities(quals);
    		}
        _visit(quals);	
        }
    
    @Override
    public void visit(final FastqRecord rec) {
    	byte quals[]= rec.getBaseQualityString().getBytes();
    	if(quals==null || quals.length==0) return;
    	SAMUtils.fastqToPhred(quals);
        _visit(quals);	
        }
   private  void _visit(final byte quals[]) {
        if(quals==null || quals.length==0) return;
        
        while(this.pos2qual.size()<quals.length)
	    	{
	    	this.pos2qual.add(new QualData());
	    	}
        
    	for(int x=0;x< quals.length;++x)
    		{
    		this.pos2qual.get(x).visit(quals[x]);
    		}
    	}
    
    public LineChart<Number, Number> build()
        {
        
    	final NumberAxis xAxis = new NumberAxis();
    	xAxis.setLabel("Position in Read");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Read Quality");

        

    	final  XYChart.Series<Number,Number> serie_mean = new XYChart.Series<Number,Number>();
    	serie_mean.setName("Mean Quality");
    	final  XYChart.Series<Number,Number> serie_min = new XYChart.Series<Number,Number>();
    	serie_min.setName("Min Quality");
    	final  XYChart.Series<Number,Number> serie_max = new XYChart.Series<Number,Number>();
    	serie_max.setName("Max Quality");
    	
        for(int x=0;x < this.pos2qual.size();++x)
        	{
        	final QualData q= this.pos2qual.get(x);
        	serie_mean.getData().add(new XYChart.Data<Number,Number>(
        			(x+1),q.get())
        			);
        	serie_min.getData().add(new XYChart.Data<Number,Number>(
        			(x+1),q.min())
        			);
        	serie_max.getData().add(new XYChart.Data<Number,Number>(
        			(x+1),q.max())
        			);
        	}
        
        
        final LineChart<Number, Number> sbc =
                new LineChart<Number, Number>(xAxis, yAxis);
        sbc.setTitle("Position/Quality");
        sbc.getData().add(serie_min);
        sbc.getData().add(serie_mean);
        sbc.getData().add(serie_max);
        sbc.setCreateSymbols(false);
        return sbc;
        }
    }
