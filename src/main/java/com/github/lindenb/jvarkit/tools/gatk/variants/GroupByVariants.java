package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import com.github.lindenb.jvarkit.gatk.Category;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * 
 *
 * 
 * 
 * 
 */
@DocumentedGATKFeature(
		summary="Reads a VCF file and creates a variant summary table",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class} )
public class GroupByVariants 
	extends AbstractGroupBy
	{    
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
    @Argument(fullName="tsv",shortName="tsv",required=false,doc="Group by Transition/Transversion")
    public boolean byTsv = false;
    @Argument(fullName="allelesize",shortName="allelesize",required=false,doc="Group by Max(allele.size)")
    public boolean byAlleleSize = false;
    @Argument(fullName="allelefrequency",shortName="af",required=false,doc="Group by Allele Frequency using the AF attribute. The lowest AF is used for multiallelic variants.")
    public boolean byAlleFrequency = false;  
    @Argument(fullName="depth",shortName="dp",required=false,doc="Group by Depth using the DP attribute")
    public boolean byDepth = false;
    @Argument(fullName="attribute",shortName="attribute",required=false,doc="Search for the presence (true/false) of an attribute in the INFO format. For example using GATK.VariantAnnotator and -comp")
    public Set<String> presence_of_attributes = new HashSet<>();
    public boolean bySingleton = false;  
    @Argument(fullName="singleton",shortName="singleton",required=false,doc="Singleton variants")

    
    private final DoubleRangeClassifier groupByAfClassifier = new DoubleRangeClassifier(
    		Arrays.asList(new DoubleRangeClassifier.Window(0.0, 0.01),new DoubleRangeClassifier.Window(0.1, 0.1))
    		, 1.0, new DecimalFormat("#.##"));
    private final IntRangeClassifier groupByDpClassifier = new IntRangeClassifier(
    		Arrays.asList(
    				new IntRangeClassifier.Window(0,1),
    				new IntRangeClassifier.Window(10,5),
    				new IntRangeClassifier.Window(50,10),
    				new IntRangeClassifier.Window(100,100))
    		, 1000,NumberFormat.getIntegerInstance());
    private final IntRangeClassifier alleleSizeClassifier = new IntRangeClassifier(
    		Arrays.asList(
    				new IntRangeClassifier.Window(0,1),
    				new IntRangeClassifier.Window(10,10),
    				new IntRangeClassifier.Window(100,100),
    				new IntRangeClassifier.Window(1000,1000)
    				)
    		, 5000,NumberFormat.getIntegerInstance());

    
    
    
	@Override
	public Map<Category,Long> map(
				final RefMetaDataTracker tracker,
				final ReferenceContext ref,
				final AlignmentContext context
				) {
		if ( tracker == null )return Collections.emptyMap();
		final Map<Category,Long> count = new HashMap<>();
		for(final VariantContext ctx: tracker.getValues(this.variants,context.getLocation()))
			{
			final List<Object> labels=new ArrayList<>();
			if(bychrom) labels.add(ctx.getContig());
			if(byID) labels.add(ctx.hasID());
			if(byType) labels.add(ctx.getType().name());
			if(byFilter) labels.add(ctx.isFiltered());
			if(minQuality>=0) {
				labels.add(ctx.hasLog10PError() && ctx.getPhredScaledQual()>=this.minQuality ?
						".":"LOWQUAL"
						);
				}
			
			if(byImpact || bybiotype)
				{
				String biotype=null;
				AnnPredictionParser.Impact impact=null;
				for(final AnnPredictionParser.AnnPrediction pred: super.annParser.getPredictions(ctx))
					{
					final AnnPredictionParser.Impact  currImpact = pred.getPutativeImpact();
					if(impact!=null && currImpact.compareTo(impact)<0) continue;
					impact=currImpact;
					biotype= pred.getTranscriptBioType();
					}
				if(byImpact) labels.add(impact==null?".":impact);
				if(bybiotype) labels.add(biotype==null?".":biotype);
				}
			if(bynalts)labels.add(ctx.getAlternateAlleles().size());
			if(byAffected || byCalled || bySingleton) 
				{
				int nc=0;
				int ng=0;
				int nsingles=0;
				for(int i=0;i< ctx.getNSamples();++i)
					{
					final Genotype g= ctx.getGenotype(i);
					if(!(!g.isCalled() || g.isNoCall() || g.isHomRef()))
						{
						ng++;
						}
					if(g.isCalled())
						{
						nc++;			
						if(!g.isHomRef())
							{
							nsingles++;
							}
						}
					}
				if(byCalled) labels.add(nc< maxSamples?nc:"GE_"+maxSamples);
				if(byAffected) labels.add(ng< maxSamples?ng:"GE_"+maxSamples);
				if(bySingleton) labels.add(nsingles==1?"SINGLETON":".");
				}
			if(byTsv) {
				if(ctx.getType()==VariantContext.Type.SNP && ctx.getAlternateAlleles().size()==1) {
					boolean b=(VariantContextUtils.isTransition(ctx));
					labels.add(b?"Transition":"Transversion");
					}
				else
					{
					labels.add(".");
					}
				
				}	
			if(byAlleleSize) {
				// see http://stackoverflow.com/questions/41678374/
				final Predicate<Allele> afilter = new Predicate<Allele>() {
					@Override
					public boolean test(final Allele a) {
						return !(a.isNoCall() || a.isSymbolic() );
					}
				};
				final OptionalInt longest = ctx.getAlleles().stream().filter(afilter).mapToInt(new ToIntFunction<Allele>() {
						public int applyAsInt(final Allele value) {return value.length();};
					}).max();
				labels.add(longest.isPresent()?alleleSizeClassifier.apply(longest.getAsInt()):"N/A");
				}
			
			if(byAlleFrequency)
				{
				final List<Object> afs= ctx.getAttributeAsList("AF");
				if(afs.isEmpty())
					{
					labels.add("NOT_AVAILABLE");
					}
				else 
					{
					Double minaf=null;
					for(final Object o:afs)
						{
						final Double af;
						if(o==null) continue;
						if(o instanceof Double) {
							af = Double.class.cast(o);
							}
						else
							{
							try
								{
								af=Double.parseDouble(String.valueOf(o));
								}
							catch(NumberFormatException err)
								{
								logger.warn("Not a number for AF :"+o);
								continue;
								}
							}
						if(af<0.0) logger.warn("AF < 0 : "+o);
						if(af>1.0) logger.warn("AF > 1.0 : "+o);
						if(minaf==null || af.compareTo(minaf)<0) 
							{
							minaf=af;
							}
						}
					labels.add(minaf==null?"NOT_FOUND":this.groupByAfClassifier.apply(minaf));
					}
				
				}
			
			if(byDepth)
				{
				final List<Object> depths= ctx.getAttributeAsList("DP");
				if(depths.size()!=1)
					{
					if(depths.size()>1) 
						{
						logger.warn("Too many data for DP :"+depths);
						}
					labels.add("NOT_AVAILABLE");
					} 
				else 
					{
					Integer dp=null;
					final Object o  = depths.get(0);
					if(o!=null && o instanceof Integer)
						{
						dp = Integer.class.cast(o);
						}
					else
						{	
						try
							{
							int i=Integer.parseInt(String.valueOf(o));
							dp=i;
							}
						catch(NumberFormatException err)
							{
							logger.warn("Not a number for DP :"+o);
							dp=null;
							}
						}
					labels.add(this.groupByDpClassifier.apply(dp));
					}
				}
			for(final String att:this.presence_of_attributes)
				{
				labels.add(ctx.hasAttribute(att));
				}
			
			final Category cat=new Category(labels);
			Long n = count.get(cat);
			count.put(cat, n==null?1L:n+1);
				
			}
		return count;
	}

   
	@Override
	protected GATKReportTable createGATKReportTable() {
		final GATKReportTable table=new GATKReportTable(
				"Variants", "Variants "+variants.getSource(),0);
		if(bychrom) table.addColumn("CONTIG");
		if(byID) table.addColumn("IN_DBSNP");
		if(byType) table.addColumn("TYPE");
		if(byFilter) table.addColumn("FILTERED");
		if(minQuality>=0) table.addColumn("QUAL_GE_"+this.minQuality);
		if(byImpact) table.addColumn("IMPACT");
		if(bybiotype) table.addColumn("BIOTYPE");
		if(bynalts) table.addColumn("N_ALT_ALLELES");
		if(byCalled) table.addColumn("CALLED_SAMPLES");
		if(byAffected) table.addColumn("AFFECTED_SAMPLES");
		if(bySingleton) table.addColumn("SINGLETON");
		if(byTsv) table.addColumn("Ts/Tv");
		if(byAlleleSize) table.addColumn("ALLELE_SIZE");
		if(byAlleFrequency) table.addColumn("ALLELE_FREQUENCY");
		if(byDepth) table.addColumn("DEPTH");
		for(final String att:this.presence_of_attributes)
			{
			 table.addColumn("ATT:"+att);
			}
		return table;
		}
	
	
	}
