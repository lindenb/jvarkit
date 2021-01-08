/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**

BEGIN_DOC

## Input

Variants in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.

VCF header must contain a pedigree ( see VCFinjectPedigree ) or a pedigree must be defined.


END_DOC
*/

@Program(name="vcfburdenfisherh",
	description="Fisher Case /Controls per Variant",
	keywords= {"vcf","burden","fisher"},
	modificationDate="20200713",
	creationDate="20160418"
	)
public class VcfBurdenFisherH
	extends OnePassVcfLauncher
	{	
	private static final Logger LOG = Logger.build(VcfBurdenFisherH.class).make();

	@Parameter(names={"-m","--min-fisher"},description="min inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double min_fisher = 0.0 ;
	@Parameter(names={"-M","--max-fisher"},description="max inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double max_fisher = 1.0 ;
	@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Don't try to calculate things why variants already FILTERed (faster)")
	private boolean ignoreFiltered=false;
	@Parameter(names={"-p","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile=null;
	@Parameter(names={"-gtf","--gtf","--gtFiltered"},description="[20180115] Ignore FILTERed **Genotype**")
	private boolean ignore_filtered_genotype=false;

	@Parameter(names={"--attribute"},description="Name of the attribute used for INFO")
	private String burdenHFisherTag = "BurdenHFisher";
	@Parameter(names={"-F1","--filter"},description="if this value is not blank, the FILTER will be set for this variant if the fisher values are out of the bounds.")
	private String softFilterName = "BurdenHFisher";
	@Parameter(names={"-F2","--ctrlgtcase"},description="Set this FILTER if the proportion of Controls carrying a ALT allele is creater than proportion of CASES. if blank, variant is discarded.")
	private String filterCtrlgtCaseRatioStr = "CTRL_CASE_RATIO";
	@Parameter(names={"--report"},description="[20190418] save report as bed file")
	private Path bedExportPath = null;
	@Parameter(names={"-Q","--qual"},description="Overwrite QUAL column with the lowest fisher value.")
	private boolean overwrite_qual=false;

	
		
	public VcfBurdenFisherH() {
	}
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter w) {
		
		final PrintWriter report;		
				
		if(this.bedExportPath==null) {
			report = new PrintWriter(new NullOuputStream());
			} 
		else
			{
			report = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(this.bedExportPath));
			}
		
		
		final VCFHeader header = r.getHeader();

		final VCFHeader h2= new VCFHeader(header);
		
	
		final Pedigree pedigree ;
		
		try {
			pedigree = new PedigreeParser().parse(this.pedigreeFile);
		} catch(final IOException error) {
			throw new RuntimeIOException(error);
		}
		
		final Set<Sample> individualSet = pedigree.getSamplesInVcfHeader(header).
				filter(S->S.isStatusSet()).
				collect(Collectors.toSet());
				
		if(individualSet.isEmpty()) throw new IllegalArgumentException("No overlapping samples between header and pedigree.");
		

		final VCFInfoHeaderLine fisherAlleleInfoHeader = new VCFInfoHeaderLine(
				this.burdenHFisherTag ,
				VCFHeaderLineCount.A,
				VCFHeaderLineType.Float,
				"Fisher Exact Test Case/Control."
				);
		final VCFFilterHeaderLine variantFilterHeader;
		if(StringUtils.isBlank(this.softFilterName)) {
			variantFilterHeader = null;
			}
		else {
			variantFilterHeader = new VCFFilterHeaderLine(
				this.softFilterName,
				"Fisher case:control vs miss|have ALT not between "+ this.min_fisher+" and "+this.max_fisher
				);
			h2.addMetaDataLine(variantFilterHeader);
			}
		
		final VCFFilterHeaderLine filterCtrlgtCaseRatio;
		
		if(StringUtils.isBlank(this.filterCtrlgtCaseRatioStr)) {
			filterCtrlgtCaseRatio = null;
			}
		else {
			filterCtrlgtCaseRatio = new VCFFilterHeaderLine(
					this.filterCtrlgtCaseRatioStr,
					"The number of CONTROLS carrying the ALT allele is creater than the number of CASES carrying the ALT allele."
					);
			h2.addMetaDataLine(filterCtrlgtCaseRatio);
			}
		
		final VCFInfoHeaderLine fisherDetailInfoHeader = new VCFInfoHeaderLine(
				this.burdenHFisherTag + "Detail",
				VCFHeaderLineCount.A,
				VCFHeaderLineType.String,
				"Fisher Exact Test Case/Control"
				);
		h2.addMetaDataLine(fisherAlleleInfoHeader);
		h2.addMetaDataLine(fisherDetailInfoHeader);
		JVarkitVersion.getInstance().addMetaData(this, h2);
		w.writeHeader(h2);
		while(r.hasNext())
			{
			final VariantContext ctx = r.next();
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			if(this.overwrite_qual) vcb.log10PError(VariantContext.NO_LOG10_PERROR);
			vcb.rmAttribute(fisherAlleleInfoHeader.getID());
			vcb.rmAttribute(fisherDetailInfoHeader.getID());
			
			final Set<String> oldFilters = new HashSet<>(ctx.getFilters());
			if(variantFilterHeader!=null) oldFilters.remove(variantFilterHeader.getID());
			if(filterCtrlgtCaseRatio!=null) oldFilters.remove(filterCtrlgtCaseRatio.getID());
			
			
			if(this.ignoreFiltered && ctx.isFiltered())
				{
				w.add(vcb.make());
				continue;
				}
						
			
			boolean set_filter_in_range = true;
			boolean set_filter_case_ctrl_ratio = true;
			boolean found_one_alt = false;
			final List<String> infoData = new ArrayList<>(ctx.getAlleles().size());
			final List<Double> fisherValues = new ArrayList<>(ctx.getAlleles().size());
			
			for(final Allele observed_alt: ctx.getAlternateAlleles()) {
				if(observed_alt.isNoCall()) {
					infoData.add(
							String.join("|",
							"ALLELE",String.valueOf(observed_alt.getDisplayString()),
							"FISHER","1.0"
							));
					fisherValues.add(1.0);
					continue;
					}
				
				/* count for fisher allele */
				int count_case_have_alt =0;
				int count_case_miss_alt = 0;
				int count_ctrl_have_alt = 0;
				int count_ctrl_miss_alt = 0;
					
				
				/* loop over persons in this pop */
				for(final Sample p: individualSet ) 	{
					/* get genotype for this individual */
					final Genotype genotype = ctx.getGenotype(p.getId());
					
					/* individual is not in vcf header */
					if(genotype==null || !genotype.isCalled() || (this.ignore_filtered_genotype && genotype.isFiltered())) {
						if(genotype==null) LOG.warn("Genotype is null for sample "+p.getId()+" not is pedigree!");
						//no information , we consider that sample was called AND HOM REF
						if(p.isAffected()) { count_case_miss_alt++; }
						else { count_ctrl_miss_alt++; }
						continue;
						}
					
					/* loop over alleles */
					final boolean genotype_contains_allele = genotype.getAlleles().stream().
							anyMatch(A->A.equals(observed_alt));
						
					/* fisher */
					if(genotype_contains_allele) {
						if(p.isAffected()) { count_case_have_alt++; ;}
						else { count_ctrl_have_alt++; }
						}
					else {
						if(p.isAffected()) { count_case_miss_alt++; }
						else { count_ctrl_miss_alt++; }
						}
				}/* end of loop over persons */
				
				
				/* fisher test for alleles */
				final FisherExactTest fisherAlt = FisherExactTest.compute(
						count_case_have_alt, count_case_miss_alt,
						count_ctrl_have_alt, count_ctrl_miss_alt
						);
				
				fisherValues.add(fisherAlt.getAsDouble());
				infoData.add(
						String.join("|",
						"ALLELE",String.valueOf(observed_alt.getDisplayString()),
						"FISHER",String.valueOf(fisherAlt.getAsDouble()),
						"CASE_HAVE_ALT",String.valueOf(count_case_have_alt),
						"CASE_MISS_ALT",String.valueOf(count_case_miss_alt),
						"CTRL_HAVE_ALT",String.valueOf(count_ctrl_have_alt),
						"CTRL_MISS_ALT",String.valueOf(count_ctrl_miss_alt)
						));
				found_one_alt = true;
				
				
				final boolean is_in_range = this.min_fisher<=fisherAlt.getAsDouble() && fisherAlt.getAsDouble() <=this.max_fisher;
				
				
				final int total_ctrls = count_ctrl_have_alt + count_ctrl_miss_alt;
				final int total_cases = count_case_have_alt + count_case_miss_alt;
				
				// check ratio case/control
				if( total_ctrls>0 && total_cases>0 && 
					(count_case_have_alt/(double)total_cases) >= (count_ctrl_have_alt/(double)total_ctrls)) {
					set_filter_case_ctrl_ratio = false;
					}
			
				if( is_in_range ) {
					set_filter_in_range = false;
				}

				
				
				
				if( this.bedExportPath != null && is_in_range ) {
					report.print(ctx.getContig());
					report.print('\t');
					report.print(ctx.getStart()-1);
					report.print('\t');
					report.print(ctx.getEnd());
					report.print('\t');
					report.print(ctx.getReference().getDisplayString());
					report.print('\t');
					report.print(observed_alt.getDisplayString());
					report.print('\t');
					report.print(fisherAlt.getAsDouble());
					report.print('\t');
					report.print(count_case_have_alt);
					report.print('\t');
					report.print(count_case_miss_alt);
					report.print('\t');
					report.print(count_ctrl_have_alt);
					report.print('\t');
					report.print(count_ctrl_miss_alt);
					report.println();
					}	
			} //end of for each ALT allele
			
			// better than nothing, otherwise i'll mess with the filters
			if(!found_one_alt) {
				w.add(vcb.make());
				continue;
			}
			
			//soft filter, skip variant
			if(
				(set_filter_in_range && variantFilterHeader==null) ||
				(set_filter_case_ctrl_ratio && filterCtrlgtCaseRatio==null)
				) continue;//skip variant
			
			vcb.attribute(fisherAlleleInfoHeader.getID(),fisherValues);
			vcb.attribute(fisherDetailInfoHeader.getID(),infoData );
			
			if(this.overwrite_qual) {
				final OptionalDouble minV = fisherValues.stream().mapToDouble(V->V.doubleValue()).min();
				if(minV.isPresent()) vcb.log10PError(Math.max(1.0E-100/* arbitrary */,minV.getAsDouble())/-10);
				}
			
			if( set_filter_in_range  && variantFilterHeader!=null) {
				vcb.filter(variantFilterHeader.getID());
				// only set this one if the filter is set above
				if( set_filter_case_ctrl_ratio  && filterCtrlgtCaseRatio!=null) {
					vcb.filter(filterCtrlgtCaseRatio.getID());
					}
				}
			
			
			w.add(vcb.make());
			}
		w.close();
		report.flush();
		report.close();
		return 0;
		}

	
	@Override
	protected int beforeVcf() {
		if(StringUtils.isBlank(this.burdenHFisherTag)) {
			LOG.error("empty  burdenHFisherAttr");
			return -1;
			}
		if(this.min_fisher>=this.max_fisher) {
			LOG.error("bad min/max.");
			return -1;
			}
		return 0;
		}
	
	
	public static void main(final String[] args)
		{
		new VcfBurdenFisherH().instanceMainWithExit(args);
		}
	}
