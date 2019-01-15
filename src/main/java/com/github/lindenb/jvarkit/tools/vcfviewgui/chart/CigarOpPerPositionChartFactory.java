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

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;

public class CigarOpPerPositionChartFactory extends ReadChartFactory
	{
    /** cigar / base / number */
    private final List<Counter<CigarOperator>> cigar2pos2count= new ArrayList<>();
    
    @Override
    public String getName() {
    	return "Cigar/Position";
    	}
    
    @Override
    public void visit(final SAMRecord rec) {
    	if(rec.getReadUnmappedFlag() || rec.getCigar()==null) return;
    	int readpos=0;
    	for(final CigarElement ce:rec.getCigar())
    		{
    		switch(ce.getOperator())
				{
				case P:break;
				case D:case N:
					{
					_visit(ce.getOperator(),readpos);
					readpos++;
					break;
					}
				default:
					for(int i=0;i< ce.getLength();++i)
						{
						_visit(ce.getOperator(),readpos);
						readpos++;
						}
					break;
					}
    			}
        }
    
    @Override
    public void visit(final FastqRecord rec) {
    	
    	}
    
   private void _visit(final CigarOperator op,int readpos) {
        while(this.cigar2pos2count.size()<=readpos)
	  		{
			this.cigar2pos2count.add(new Counter<>());
	  		}
        this.cigar2pos2count.get(readpos).incr(op);
    	}
   
    @Override
    public StackedBarChart<String, Number> build()
        {        
    	final CategoryAxis xAxis = new CategoryAxis();
    	xAxis.setLabel("Position in Read");
        final NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("Count");

        
        final List<XYChart.Series<String, Number>> base2count=new ArrayList<>();
        for(final CigarOperator cigarop:CigarOperator.values())
        	{
        	if(cigarop==CigarOperator.P) continue;
        	final XYChart.Series<String, Number> serie= new XYChart.Series<String, Number>();
        	serie.setName(cigarop.name());
        	base2count.add(serie);
        	
	        for(int i=0;i<  this.cigar2pos2count.size();++i)
	        	{
	        	serie.getData().add(new XYChart.Data<String,Number>(
	        			String.valueOf(i+1),
	        			this.cigar2pos2count.get(i).count(cigarop))
	        			);
	        	}
	        }
        
        final StackedBarChart<String, Number> sbc =
                new StackedBarChart<String, Number>(xAxis, yAxis);
        sbc.setTitle(getName());
        sbc.getData().addAll(base2count);
        sbc.setCategoryGap(0.2);
        
        return sbc;
        }
    }
