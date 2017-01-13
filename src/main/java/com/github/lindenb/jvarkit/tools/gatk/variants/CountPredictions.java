package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Test Documentation
 *
 * 
 * 
 * 
 */
@DocumentedGATKFeature(
		summary="Count Predictions",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class} )
public class CountPredictions  extends RodWalker<Map<CountPredictions.Category,Long>, Map<CountPredictions.Category,Long>> {
	private enum IMPACT { HIGH, MODERATE, MODIFIER,LOW};
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    @Output(doc="File to which result should be written")
    protected PrintStream out = System.out;
    
    @Argument(fullName="minQuality",shortName="mq",required=false,doc="Group by Quality. Set the treshold for Minimum Quality")
    public double minQuality = -1;
    @Argument(fullName="chrom",shortName="chrom",required=false,doc="Group by Chromosome/Contig")
    public boolean bychrom = false;
    @Argument(fullName="ID",shortName="ID",required=false,doc="Group by ID")
    public boolean byID = false;
    @Argument(fullName="variantType",shortName="variantType",required=false,doc="Group by VariantType")
    public boolean byType = false;
    @Argument(fullName="filter",shortName="filter",required=false,doc="Group by FILTER")
    public boolean byFilter = false;
    @Argument(fullName="impact",shortName="impact",required=false,doc="Group by ANN/IMPACT")
    public boolean byImpact = false;
    @Argument(fullName="biotype",shortName="biotype",required=false,doc="Group by ANN/biotype")
    public boolean bybiotype = false;
    @Argument(fullName="nalts",shortName="nalts",required=false,doc="Group by number of ALTS")
    public boolean bynalts = false;
    @Argument(fullName="affected",shortName="affected",required=false,doc="Group by number of Samples called and not HOMREF")
    public boolean byAffected = false;
    @Argument(fullName="called",shortName="called",required=false,doc="Group by number of Samples called")
    public boolean byCalled = false;
    @Argument(fullName="maxSamples",shortName="maxSamples",required=false,doc="if the number of samples affected is greater than --maxSamples use the label \"GT_MAX_SAMPLES\"")
    public int maxSamples = Integer.MAX_VALUE;
   
    private final Pattern pipeRegex=Pattern.compile("[\\|]");
    
    static class Category
    	implements Comparable<Category>
    	{
    	private final List<String> labels ;
    	Category(final List<String> labels) {
    		this.labels=new ArrayList<>(labels);
    		}
    	@Override
    	public int hashCode() {
    		return labels.hashCode();
    		}
    	@Override
    	public int compareTo(final Category o) {
    		for(int i=0;i< labels.size();++i)
    			{
    			int d = labels.get(i).compareTo(o.labels.get(i));
    			if(d!=0) return d;
    			}
    		return 0;
    		}
    	@Override
    	public boolean equals(Object o) {
    		if(o==this) return true;
    		if(o==null || !(o instanceof Category)) return false;
    		return compareTo(Category.class.cast(o))==0;
    		}
    	}
    private int ann_impact_column=-1;/* eg: MODIFIER / LOW */
    private int ann_transcript_biotype_column=-1;/* eg: transcript /intergenic_region */

    @Override
    public void initialize() {
    	
    	if(byImpact || bybiotype) {
        	final VCFHeader vcfHeader  = GATKVCFUtils.getVCFHeadersFromRods(getToolkit()).get(variants.getName());
	    	final VCFInfoHeaderLine annInfo = vcfHeader.getInfoHeaderLine("ANN");
	    	if(annInfo==null)
	    		{
	    		logger.warn("NO ANN in "+variants.getSource());
	    		}
	    	else
	    		{
	    		int q0=annInfo.getDescription().indexOf('\'');
	    		int q1=annInfo.getDescription().lastIndexOf('\'');
	    		if(q0==-1 || q1<=q0)
	    			{
	    			logger.warn("Cannot parse "+annInfo.getDescription());
	    			}
	    		else
		    		{
		    		final String fields[]=pipeRegex.split(annInfo.getDescription().substring(q0+1, q1));
		    		for(int c=0;c<fields.length;++c)
		    	    	{
		    	    	final String column=fields[c].trim();
		    	    	if(column.equals("Annotation_Impact"))
				    		{
		    	    		ann_impact_column=c;
				    		}
		    	    	else if(column.equals("Transcript_BioType")) 
		    	    		{
		    	    		ann_transcript_biotype_column=c;
		    	    		}
		    	    	}
		    		}
	    		}
	    	}
    	
    	super.initialize();
    }
    
	@Override
	public Map<CountPredictions.Category,Long> map(final RefMetaDataTracker tracker,final ReferenceContext ref, final AlignmentContext context) {
		if ( tracker == null )return Collections.emptyMap();
		final Map<CountPredictions.Category,Long> count = new HashMap<>();
		for(final VariantContext ctx: tracker.getValues(this.variants,context.getLocation()))
			{
			final List<String> labels=new ArrayList<>();
			if(bychrom) labels.add(ctx.getContig());
			if(byID) labels.add(ctx.hasID()?"Y":".");
			if(byType) labels.add(ctx.getType().name());
			if(byFilter) labels.add(ctx.isFiltered()?"F":".");
			if(minQuality>=0) {
				labels.add(ctx.hasLog10PError() && ctx.getPhredScaledQual()>=this.minQuality ?
						".":"LOWQUAL"
						);
				}
			
			if(byImpact || bybiotype)
				{
				String biotype=null;
				IMPACT impact=null;
				final List<Object> anns = ctx.getAttributeAsList("ANN");
				for(final Object anno:anns) {
					final String tokens[]=this.pipeRegex.split(anno.toString());
					if(this.ann_impact_column==-1 ||
							this.ann_impact_column >= tokens.length ||
							tokens[this.ann_impact_column].isEmpty()
							) continue;
					IMPACT currImpact = IMPACT.valueOf(tokens[this.ann_impact_column]);
					if(impact!=null && currImpact.compareTo(impact)<0) continue;
					impact=currImpact;
					biotype=null;
					if(this.ann_transcript_biotype_column==-1 ||
							this.ann_transcript_biotype_column >= tokens.length ||
							tokens[this.ann_transcript_biotype_column].isEmpty()
							) continue;
					biotype=tokens[ann_transcript_biotype_column];
					}
				if(byImpact) labels.add(impact==null?".":impact.name());
				if(bybiotype) labels.add(biotype==null?".":biotype);
				}
			if(bynalts)labels.add(String.valueOf(ctx.getAlternateAlleles().size()));
			if(byAffected || byCalled) 
				{
				int nc=0;
				int ng=0;
				for(int i=0;i< ctx.getNSamples();++i)
					{
					final Genotype g= ctx.getGenotype(i);
					if(!(g.isNoCall() || g.isHomRef()))
						{
						ng++;
						}
					if(g.isCalled())
						{
						nc++;
						}
					}
				if(byCalled) labels.add(nc< maxSamples?String.valueOf(nc):"GE_"+maxSamples);
				if(byAffected) labels.add(ng< maxSamples?String.valueOf(ng):"GE_"+maxSamples);
				}
			
			final Category cat=new Category(labels);
			Long n=count.get(cat);
			count.put(cat, n==null?1L:n+1);
				
			}
		return count;
	}

	@Override
	public Map<CountPredictions.Category,Long> reduce(Map<CountPredictions.Category,Long> value, Map<CountPredictions.Category,Long> sum) {
		final Map<CountPredictions.Category,Long> newmap = new HashMap<>(sum);
		for(Category cat:value.keySet()) {
			Long sv = sum.get(cat);
			Long vv = value.get(cat);
			newmap.put(cat, sv==null?vv:sv+vv);
		}
		return newmap;
	}

	@Override
	public Map<CountPredictions.Category,Long> reduceInit() {
		return Collections.emptyMap();
	}
   

	@Override
	public void onTraversalDone(final Map<CountPredictions.Category,Long> counts) {
		GATKReportTable table=new GATKReportTable(
				"Variants", "Variants "+variants.getSource(),0);
		if(bychrom) table.addColumn("CONTIG");
		if(byID) table.addColumn("IN_DBSNP");
		if(byType) table.addColumn("TYPE");
		if(byFilter) table.addColumn("FILTER");
		if(minQuality>=0) table.addColumn("QUAL_GE_"+this.minQuality);
		if(byImpact) table.addColumn("IMPACT");
		if(bybiotype) table.addColumn("BIOTYPE");
		if(bynalts) table.addColumn("N_ALT_ALLELES");
		if(byCalled) table.addColumn("CALLED_SAMPLES");
		if(byAffected) table.addColumn("AFFECTED_SAMPLES");
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
		GATKReport report = new GATKReport();
		report.addTable(table);
		report.print(this.out);
		out.flush();
		
		logger.info("TraversalDone");

		}
	
	}
