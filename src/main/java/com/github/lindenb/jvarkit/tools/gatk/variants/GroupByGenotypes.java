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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

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

/**
 * Test Documentation
 *
 * <p>Hello world</p>
 * 
 * 
 */
@DocumentedGATKFeature(
		summary="Reads a VCF file and creates a genotype summary table",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class}
		)
public class GroupByGenotypes  extends AbstractGroupBy {
    
    @Argument(fullName="minGenotypeQuality",shortName="mgq",required=false,doc="Minimum genotype quality to put a assign a category")
    public int minGenotypeQuality = 0;
    @Argument(fullName="chrom",shortName="chrom",required=false,doc="Group by Chromosome/Contig")
    public boolean bychrom = false;
    @Argument(fullName="ID",shortName="ID",required=false,doc="Group by ID")
    public boolean byID = false;
    @Argument(fullName="variantType",shortName="variantType",required=false,doc="Group by VariantType")
    public boolean byType = false;
    @Argument(fullName="genotypeType",shortName="genotypeType",required=false,doc="Group by GenotypeType")
    public boolean byGenotypeType = false;

    @Argument(fullName="filter",shortName="filter",required=false,doc="Group by FILTER")
    public boolean byFilter = false;
    @Argument(fullName="gfilter",shortName="gfilter",required=false,doc="Group by GENOTYPE FILTER")
    public boolean byGFilter = false;

    @Argument(fullName="impact",shortName="impact",required=false,doc="Group by ANN/IMPACT")
    public boolean byImpact = false;
    @Argument(fullName="onlysingletons",shortName="onlysingletons",required=false,doc="only consider singletons (one sample affected/variant)")
    public boolean onlysingletons = false;

    
    @Override
    public Map<Category, Long> map(
    		final RefMetaDataTracker tracker, 
    		final ReferenceContext refctx,
    		final AlignmentContext context
    		) {
    	if ( tracker == null )return Collections.emptyMap();
        final Map<Category, Long> counts = new HashMap<>();

		for(final VariantContext ctx: tracker.getValues(this.variants,context.getLocation()))
			{
			int index_singleton=-1;
			if(onlysingletons) {
				for(int i=0;i< ctx.getNSamples();++i)
					{
					final Genotype g = ctx.getGenotype(i);
					if(g==null || !g.isCalled() || g.isNoCall() || g.isHomRef()) continue;
					if(index_singleton!=-1) {
						//not anymore a singleton
						index_singleton=-1;
						break;
						}
					index_singleton = i;
					}
				}
			for(int i=0;i< ctx.getNSamples();++i)
				{
				
				if(onlysingletons && index_singleton!=i) {
					continue;
					}
				final Genotype genotype=ctx.getGenotype(i);
				final List<Object> labels=new ArrayList<>();
				labels.add(genotype.getSampleName());
				if(bychrom) labels.add(ctx.getContig());
				if(byID) labels.add(ctx.hasID());
				if(byType) labels.add(ctx.getType().name());
				if(byGenotypeType) labels.add(genotype.getType());
				if(byFilter) labels.add(ctx.isFiltered());
				if(byGFilter) labels.add(genotype.isFiltered());
				
				if(minGenotypeQuality>=0)
					{
					labels.add(genotype.hasGQ() && genotype.getGQ()>=this.minGenotypeQuality ?
						".":"LOWQUAL"
						);
					}
				
				if(byImpact)
					{
					AnnPredictionParser.Impact impact=null;
					for(final AnnPredictionParser.AnnPrediction pred: super.annParser.getPredictions(ctx)) {
						// see http://stackoverflow.com/questions/41678374/
						final Predicate<Allele> afilter = new Predicate<Allele>() {							
							@Override
							public boolean test(final Allele A) {
								return A.getDisplayString().equals(pred.getAllele());
							}
						};
						
						if(genotype.getAlleles().
								stream().
								filter(afilter).
								findAny().isPresent()==false) continue;
						
						final AnnPredictionParser.Impact currImpact = pred.getPutativeImpact();
						if(impact!=null && currImpact.compareTo(impact)<0) continue;
						impact=currImpact;
						}
					if(byImpact) labels.add(impact==null?".":impact.name());
					}
				
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
		if(bychrom) table.addColumn("CHROM");
		if(byID)table.addColumn("IN_DBSNP");
		if(byType) table.addColumn("VARIANT_TYPE");
		if(byGenotypeType) table.addColumn("GENOTYPE_TYPE");
		if(byFilter) table.addColumn("VARIANT_FILTERED");
		if(byGFilter) table.addColumn("GENOTYPE_FILTERED");	
		if(minGenotypeQuality>=0) table.addColumn("GENOTYPE_QUAL_GE_"+this.minGenotypeQuality);
		if(byImpact) table.addColumn("IMPACT");
		
	return table;	
	}

	
	}
