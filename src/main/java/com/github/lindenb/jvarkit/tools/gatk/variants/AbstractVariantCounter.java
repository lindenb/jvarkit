package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.PrintStream;
import java.util.ArrayList;
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
 * Test Documentation
 *
 * 
 * 
 * 
 */
public abstract class AbstractVariantCounter 
	extends RodWalker<Map<AbstractVariantCounter.Category,Long>, Map<AbstractVariantCounter.Category,Long>> 
	implements  org.broadinstitute.gatk.engine.walkers.TreeReducible< Map<AbstractVariantCounter.Category,Long> >
	{
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    @Output(doc="File to which result should be written")
    public PrintStream out = System.out;
        
    protected AnnPredictionParser annParser= new PredictionParserFactory().buildAnnPredictionParser();
    
    static class Category
    	{
    	private final List<Object> labels ;
    	Category(final List<Object> labels) {
    		this.labels=new ArrayList<>(labels);
    		}
    	@Override
    	public int hashCode() {
    		return labels.hashCode();
    		}
    	@Override
    	public boolean equals(Object o) {
    		if(o==this) return true;
    		if(o==null || !(o instanceof Category)) return false;
    		return this.labels.equals(Category.class.cast(o).labels);
    		}
    	}
    @Override
    public void initialize() {
    	
        	final VCFHeader vcfHeader  = GATKVCFUtils.getVCFHeadersFromRods(getToolkit()).get(variants.getName());
	    	this.annParser = new PredictionParserFactory().header(vcfHeader).buildAnnPredictionParser();
	    	
    	
    	super.initialize();
    }
    
	
	@Override
	public Map<AbstractVariantCounter.Category,Long> treeReduce(Map<AbstractVariantCounter.Category,Long> value, Map<AbstractVariantCounter.Category,Long> sum) {
		return this.reduce(value,sum);
		}

	@Override
	public Map<AbstractVariantCounter.Category,Long> reduce(Map<AbstractVariantCounter.Category,Long> value, Map<AbstractVariantCounter.Category,Long> sum) {
		final Map<AbstractVariantCounter.Category,Long> newmap = new HashMap<>(sum);
		for(final Category cat:value.keySet()) {
			final Long sv = sum.get(cat);
			final Long vv = value.get(cat);
			newmap.put(cat, sv==null?vv:sv+vv);
		}
		return newmap;
	}

	@Override
	public Map<AbstractVariantCounter.Category,Long> reduceInit() {
		return Collections.emptyMap();
	}
	protected abstract  GATKReportTable createGATKReportTable();

	@Override
	public void onTraversalDone(final Map<AbstractVariantCounter.Category,Long> counts) {
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
		out.flush();
		
		logger.info("TraversalDone");

		}
	
	}
