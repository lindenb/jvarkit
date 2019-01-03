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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import com.github.lindenb.jvarkit.gatk.Category;
import com.github.lindenb.jvarkit.gatk.GatkReportWriter;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
 * AbstractGroupBy
 * 
 */
public abstract class AbstractGroupBy
	extends RodWalker<Map<Category,Long>, Map<Category,Long>> 
	implements  org.broadinstitute.gatk.engine.walkers.TreeReducible< Map<Category,Long> >
	{
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    @Output(doc="File to which result should be written")
    public PrintStream out = System.out;
    @Argument(shortName="outputTableFormat",fullName="outputTableFormat", doc="Format for the output table", required=false)
    protected GatkReportWriter.Format outputTableFormat =  GatkReportWriter.Format.DEFAULT;

    
    protected AnnPredictionParser annParser= new AnnPredictionParserFactory().get();
    
    /** utility class to get a window/range result for doubles values */
    protected static class DoubleRangeClassifier implements Function<Object, String>
    	{
    	static class Window
    		{
    		private final double min;
    		private final double window_size;
    		Window(
        			final double min,
        			final double window_size
        			)
        		{
        		this.min=min;
        		this.window_size=window_size;
        		}
    		@Override
    		public String toString() {
    			return "(min:"+min+" size:"+window_size+")";
    			}
    		}
    	
    	private final List<Window> windows;
    	private final double max;
    	private final DecimalFormat decimalFormat;
    	DoubleRangeClassifier(
    			List<Window> windows,
    			final double max,
    			final DecimalFormat decimalFormat
    			)
    		{
    		List<Window> newins =new ArrayList<>(windows);
    		newins.sort(new Comparator<Window>() {
	    			@Override
	    			public int compare(final Window o1,final Window o2) {
	    				if(o1.min<o2.min) return -1;
	    				if(o1.min>o2.min) return 1;
	    				throw new IllegalStateException("got to windows with same min value");
	    			}
	    		});
    		this.windows=Collections.unmodifiableList(newins);
    		if(this.windows.isEmpty()) throw new IllegalArgumentException("No window defined");
    		this.max=max;
    		for(final Window w:this.windows) if(this.max<=w.min)  throw new IllegalArgumentException("Windows min>=max "+w);
    		this.decimalFormat = decimalFormat;
    		}
    	@Override
    	public String apply(final Object t) {
    		if(t==null) return "N/A";
    		else if(t instanceof Double) 
				{
    			final double v =Double.class.cast(t);
				if(v<this.windows.get(0).min) return "LT_"+this.decimalFormat.format(this.windows.get(0).min);
				if(v>=this.max) return "GE_"+this.decimalFormat.format(this.max);

    			for(int i=0;i< windows.size();++i)
    				{
    				//double M = this.max;
    				if(i+1<windows.size())
    					{
    					if( v >=windows.get(i+1).min) continue;
    					//M = windows.get(i+1).min;
    					}
    				final Window win = windows.get(i);
    				if(v<win.min) continue;
    				double x= win.min;
    				for(;;)
    					{
    					if(x<= v && v< (x+win.window_size)) {
							return "[" +
								this.decimalFormat.format(x) +
								"-" +
								this.decimalFormat.format(x+win.window_size) +
								"["
								;
    						}
    					x+=win.window_size;
    					}    				
    				}
    			throw new IllegalStateException("should never happen "+v + " "+this.windows);
				}
    		else
    			{
    			try {
					return this.apply(Double.parseDouble(t.toString()));
				} catch (NumberFormatException e) {
					return "NaN";
					}
    			}
    		}
    	}
    
    /** utility class to get a window/range result for integer values */
    protected static class IntRangeClassifier implements Function<Object, String>
    {
    	static class Window
    	{
    		private final int min;
    		private final int window_size;
    		Window(
    				final int min,
    				final int window_size
    				)
    		{
    			this.min=min;
    			this.window_size=window_size;
    		}
    		@Override
    		public String toString() {
    			return "(min:"+min+" size:"+window_size+")";
    			}
    	}

    	private final List<Window> windows;
    	private final int max;
    	private final NumberFormat decimalFormat;
    	IntRangeClassifier(
    			List<Window> windows,
    			final int max,
    			final NumberFormat decimalFormat
    			)
    	{
    		List<Window> newins =new ArrayList<>(windows);
    		newins.sort(new Comparator<Window>() {
    			@Override
    			public int compare(final Window o1,final Window o2) {
    				if(o1.min<o2.min) return -1;
    				if(o1.min>o2.min) return 1;
    				throw new IllegalStateException("got to windows with same min value");
    			}
    		});
    		this.windows=Collections.unmodifiableList(newins);
    		if(this.windows.isEmpty()) throw new IllegalArgumentException("No window defined");
    		this.max=max;
    		for(final Window w:this.windows) if(this.max<=w.min)  throw new IllegalArgumentException("Windows min>=max "+w);
    		this.decimalFormat = decimalFormat;
    	}
    	@Override
    	public String apply(final Object t) {
    		if(t==null) return "N/A";
    		else if(t instanceof Integer) 
    			{
    			final int v =Integer.class.cast(t);
    			if(v<this.windows.get(0).min) return "LT_"+this.decimalFormat.format(this.windows.get(0).min);
    			if(v>=this.max) return "GE_"+this.decimalFormat.format(this.max);
    			for(int i=0;i< windows.size();++i)
    			{
    				//int M = this.max;
    				
    				if(i+1<windows.size())
    				{
    					if( v >=windows.get(i+1).min) continue;
    					//M = windows.get(i+1).min;
    				}
    				final Window win = windows.get(i);
    				if(v<win.min) continue;
    				int x= win.min;
    				
    				for(;;)
    					{
    					//if(v==944) System.err.println("x"+x+" v="+v+" "+ win);
    					if(x<= v && v< (x+win.window_size)/* && (x+win.window_size)<M*/) 
    						{	
    						if(win.window_size==1)
    							{
    							return this.decimalFormat.format(x);
    							}
    						else
    							{
    							return "[" +
									this.decimalFormat.format(x) +
									"-" +
									this.decimalFormat.format(x+win.window_size) +
									"["
									;
    							}
    						}
	    				x+=win.window_size;
	    				}    				
	    			}
    			throw new IllegalStateException("should never happen v="+v+" windows="+this.windows);
    		}
    		else
    		{
    			try {
    				return this.apply(Integer.parseInt(t.toString()));
    			} catch (NumberFormatException e) {
    				return "NaN";
    			}
    		}
    	}
    }
    
    
    @Override
    public void initialize() {
    	final VCFHeader vcfHeader  = GATKVCFUtils.getVCFHeadersFromRods(getToolkit()).get(variants.getName());
    	this.annParser = new AnnPredictionParserFactory(vcfHeader).get();
    	super.initialize();
    	}
    
	
	@Override
	public Map<Category,Long> treeReduce(
			final Map<Category,Long> value,
			final Map<Category,Long> sum) {
		return this.reduce(value,sum);
		}

	@Override
	public Map<Category,Long> reduce(
		final Map<Category,Long> value,
		final Map<Category,Long> sum) {
		final Map<Category,Long> newmap = new HashMap<>(sum);
		for(final Category cat:value.keySet()) {
			final Long sv = sum.get(cat);
			final Long vv = value.get(cat);
			newmap.put(cat, sv==null?vv:sv+vv);
		}
		return newmap;
	}

	@Override
	public Map<Category,Long> reduceInit() {
		return Collections.emptyMap();
	}
	protected abstract  GATKReportTable createGATKReportTable();

	@Override
	public void onTraversalDone(final Map<Category,Long> counts) {
		final GATKReportTable table=createGATKReportTable();
		table.addColumn("COUNT");
		
		int nRows=0;
		for(final Category cat: counts.keySet())
			{
			for(int x=0;x<cat.size();++x)
				{
				table.set(nRows, x, cat.get(x));
				}
			table.set(nRows, cat.size(), counts.get(cat));
			++nRows;
			}
		final GatkReportWriter reportWriter = GatkReportWriter.createWriter(this.outputTableFormat);
		final GATKReport report = new GATKReport();
		report.addTable(table);
		reportWriter.print(report, this.out);
		this.out.flush();
		
		logger.info("TraversalDone");
		}
	
	
	}
