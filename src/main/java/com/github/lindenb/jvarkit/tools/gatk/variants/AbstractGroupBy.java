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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.PrintStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.PredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
 * AbstractGroupBy
 * 
 */
abstract class AbstractGroupBy 
	extends RodWalker<Map<AbstractGroupBy.Category,Long>, Map<AbstractGroupBy.Category,Long>> 
	implements  org.broadinstitute.gatk.engine.walkers.TreeReducible< Map<AbstractGroupBy.Category,Long> >
	{
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    @Output(doc="File to which result should be written")
    public PrintStream out = System.out;
        
    protected AnnPredictionParser annParser= new PredictionParserFactory().buildAnnPredictionParser();
    
    static public class Category
    	{
    	private final int _hash;
    	private final List<Object> labels ;
    	Category(final List<Object> labels) {
    		this.labels= Collections.unmodifiableList(labels);
    		this._hash = this.labels.hashCode();
    		}
    	@Override
    	public int hashCode() {
    		return this._hash;
    		}
    	@Override
    	public boolean equals(Object o) {
    		if(o==this) return true;
    		if(o==null || !(o instanceof Category)) return false;
    		return this.labels.equals(Category.class.cast(o).labels);
    		}
    	}
    
    protected AbstractGroupBy()
    {
    	
    }
    
    @Override
    public void initialize() {
    	
        	final VCFHeader vcfHeader  = GATKVCFUtils.getVCFHeadersFromRods(getToolkit()).get(variants.getName());
	    	this.annParser = new PredictionParserFactory().header(vcfHeader).buildAnnPredictionParser();
	    	
    	
    	super.initialize();
    }
    
	
	@Override
	public Map<AbstractGroupBy.Category,Long> treeReduce(Map<AbstractGroupBy.Category,Long> value, Map<AbstractGroupBy.Category,Long> sum) {
		return this.reduce(value,sum);
		}

	@Override
	public Map<AbstractGroupBy.Category,Long> reduce(Map<AbstractGroupBy.Category,Long> value, Map<AbstractGroupBy.Category,Long> sum) {
		final Map<AbstractGroupBy.Category,Long> newmap = new HashMap<>(sum);
		for(final Category cat:value.keySet()) {
			final Long sv = sum.get(cat);
			final Long vv = value.get(cat);
			newmap.put(cat, sv==null?vv:sv+vv);
		}
		return newmap;
	}

	@Override
	public Map<AbstractGroupBy.Category,Long> reduceInit() {
		return Collections.emptyMap();
	}
	protected abstract  GATKReportTable createGATKReportTable();

	@Override
	public void onTraversalDone(final Map<AbstractGroupBy.Category,Long> counts) {
		final GATKReportTable table=createGATKReportTable();
		table.addColumn("COUNT");
		
		int nRows=0;
		for(final Category cat: counts.keySet())
			{
			for(int x=0;x<cat.labels.size();++x)
				{
				table.set(nRows, x, cat.labels.get(x));
				}
			table.set(nRows, cat.labels.size(), counts.get(cat));
			++nRows;
			}
		final GATKReport report = new GATKReport();
		report.addTable(table);
		report.print(this.out);
		this.out.flush();
		
		logger.info("TraversalDone");
		}
	
	}
