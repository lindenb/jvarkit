package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Test Documentation
 *
 * <p>Hello world</p>
 * 
 * 
 */
@DocumentedGATKFeature(
		summary="Count Variants",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class} )
public class CountSingletons  extends RodWalker<Integer, Integer> {
	private enum IMPACT { HIGH, MODERATE, MODIFIER,LOW};
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    
    @Argument(fullName="minGenotypeQuality",shortName="mgq",required=false,doc="Minimum genotype quality to put a assign a category")
    public int minGenotypeQuality = 0;

    
    private final Pattern pipeRegex=Pattern.compile("[\\|]");
    
    private static class Category
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
    private int ann_allele_column=-1;/* eg: A ... */
    private int ann_annotation_column=-1;/* eg: downstream_gene_variant ... */
    private int ann_impact_column=-1;/* eg: MODIFIER / LOW */
    private int ann_transcript_biotype_column=-1;/* eg: transcript /intergenic_region */

    private Map<Category, Long> counts = new HashMap<>();
    @Override
    public void initialize() {
    	final VCFHeader vcfHeader  = GATKVCFUtils.getVCFHeadersFromRods(getToolkit()).get(variants.getName());
    	if(vcfHeader.getNGenotypeSamples()==0)
		{
			throw new UserException("No sample in "+variants.getSource());
		}
    	final VCFInfoHeaderLine annInfo = vcfHeader.getInfoHeaderLine("ANN");
    	if(annInfo==null)
    		{
    		super.logger.warn("NO ANN in "+variants.getSource());
    		}
    	else
    		{
    		int q0=annInfo.getDescription().indexOf('\'');
    		int q1=annInfo.getDescription().lastIndexOf('\'');
    		if(q0==-1 || q1<=q0)
    			{
    			super.logger.warn("Cannot parse "+annInfo.getDescription());
    			}
    		else
	    		{
	    		final String fields[]=pipeRegex.split(annInfo.getDescription().substring(q0+1, q1));
	    		for(int c=0;c<fields.length;++c)
	    	    	{
	    	    	final String column=fields[c].trim();
	    	    	if(column.equalsIgnoreCase("Allele"))
	    	    		{
	    	    		ann_allele_column=c;
	    	    		}
	    	    	else if(column.equals("Annotation"))
			    		{
	    	    		ann_annotation_column=c;
			    		}
	    	    	else if(column.equals("Annotation_Impact"))
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
    	
    	this.counts.clear();
    	super.initialize();
    }
    
	@Override
	public Integer map(final RefMetaDataTracker tracker,final ReferenceContext ref, final AlignmentContext context) {
		if ( tracker == null )return 0;
		for(final VariantContext ctx: tracker.getValues(this.variants,context.getLocation()))
			{
			Genotype singletonGenotype=null;
			for(int i=0;i< ctx.getNSamples();++i)
				{
				final Genotype g = ctx.getGenotype(i);
				if(g==null || g.isNoCall() || g.isHomRef()) continue;
				if(singletonGenotype!=null) {
					//not anymore a singleton
					singletonGenotype=null;
					break;
					}
				singletonGenotype = g;
				}
			if(singletonGenotype!=null) {
				IMPACT impact=null;
				final List<String> labels=new ArrayList<>();
				labels.add(singletonGenotype.getSampleName());
				labels.add(ctx.hasID()?"Y":".");
				labels.add(ctx.getType().name());
				labels.add(ctx.isFiltered() || singletonGenotype.isFiltered() ?"BAD":".");
				labels.add(singletonGenotype.hasGQ() && singletonGenotype.getGQ()>=this.minGenotypeQuality ?
					".":"LOWQUAL"
					);
				
				
				final List<Object> anns = ctx.getAttributeAsList("ANN");
				for(final Object anno:anns) {
					final String tokens[]=this.pipeRegex.split(anno.toString());
					if(this.ann_allele_column==-1 ||
						this.ann_allele_column >= tokens.length ||
						tokens[this.ann_allele_column].isEmpty()
						) continue;
					final Allele alt=Allele.create(tokens[this.ann_allele_column],false);
					if(!singletonGenotype.getAlleles().contains(alt)) continue;
					if(this.ann_impact_column==-1 ||
							this.ann_impact_column >= tokens.length ||
							tokens[this.ann_impact_column].isEmpty()
							) continue;
					IMPACT currImpact = IMPACT.valueOf(tokens[this.ann_impact_column]);
					if(impact!=null && currImpact.compareTo(impact)<0) continue;
					impact=currImpact;
					}
				labels.add(impact==null?".":impact.name());
				
				final Category cat=new Category(labels);
				Long n=counts.get(cat);
				counts.put(cat, n==null?1L:n+1);
				}
			}
		return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}
   

	@Override
	public void onTraversalDone(Integer result) {
		GATKReportTable table=new GATKReportTable(
				"Singletons", "Singletons in "+variants.getSource(),7);
		table.addColumn("SAMPLE_NAME");
		table.addColumn("IN_DBSNP");
		table.addColumn("TYPE");
		table.addColumn("FILTERED");
		table.addColumn("IMPACT");
		table.addColumn("GENOTYPE_QUAL_GE_"+this.minGenotypeQuality);
		table.addColumn("COUNT");
		int nRows=0;
		for(final Category cat:this.counts.keySet())
			{
			table.set(nRows, 0, cat.labels.get(0));
			table.set(nRows, 1, cat.labels.get(1));
			table.set(nRows, 2, cat.labels.get(2));
			table.set(nRows, 3, cat.labels.get(3));
			table.set(nRows, 4, cat.labels.get(4));
			table.set(nRows, 5, cat.labels.get(5));
			table.set(nRows, 6, counts.get(cat));
			++nRows;
			}
		GATKReport report = new GATKReport();
		report.addTable(table);
		report.print(System.out);
		
		logger.info("TraversalDone");
		}
	
	}
