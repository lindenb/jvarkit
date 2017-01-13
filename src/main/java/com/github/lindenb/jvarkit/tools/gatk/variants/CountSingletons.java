package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

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
public class CountSingletons  extends AbstractVariantCounter {
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;
    
    @Argument(fullName="minGenotypeQuality",shortName="mgq",required=false,doc="Minimum genotype quality to put a assign a category")
    public int minGenotypeQuality = 0;

	
    @Override
    public Map<com.github.lindenb.jvarkit.tools.gatk.variants.AbstractVariantCounter.Category, Long> map(
    		RefMetaDataTracker tracker, ReferenceContext refctx, AlignmentContext context) {
    	if ( tracker == null )return Collections.emptyMap();
        final Map<Category, Long> counts = new HashMap<>();

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
				AnnPredictionParser.Impact impact=null;
				final List<Object> labels=new ArrayList<>();
				labels.add(singletonGenotype.getSampleName());
				labels.add(ctx.hasID()?"Y":".");
				labels.add(ctx.getType().name());
				labels.add(ctx.isFiltered() || singletonGenotype.isFiltered() ?"BAD":".");
				labels.add(singletonGenotype.hasGQ() && singletonGenotype.getGQ()>=this.minGenotypeQuality ?
					".":"LOWQUAL"
					);
				
				
				for(final AnnPredictionParser.AnnPrediction pred: super.annParser.getPredictions(ctx)) {
					
					if(singletonGenotype.getAlleles().
							stream().
							filter(A->A.getDisplayString().equals(pred.getAllele())).
							findAny().isPresent()==false) continue;
					
					AnnPredictionParser.Impact currImpact = pred.getPutativeImpact();
					if(impact!=null && currImpact.compareTo(impact)<0) continue;
					impact=currImpact;
					}
				labels.add(impact==null?".":impact.name());
				
				final Category cat=new Category(labels);
				Long n=counts.get(cat);
				counts.put(cat, n==null?1L:n+1);
				}
			}
		return counts;
	}

   
	@Override
	protected GATKReportTable createGATKReportTable() {
		final GATKReportTable table=new GATKReportTable(
				"Singletons", "Singletons in "+variants.getSource(),7);
		table.addColumn("SAMPLE_NAME");
		table.addColumn("IN_DBSNP");
		table.addColumn("TYPE");
		table.addColumn("FILTERED");
		table.addColumn("IMPACT");
		table.addColumn("GENOTYPE_QUAL_GE_"+this.minGenotypeQuality);
	return table;	
	}

	
	}
