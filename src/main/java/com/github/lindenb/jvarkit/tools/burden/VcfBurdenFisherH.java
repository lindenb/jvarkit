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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

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
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.tools.lumpysv.LumpyConstants;
import com.github.lindenb.jvarkit.util.Pedigree;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**

BEGIN_DOC

## Input

Variants in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.

VCF header must contain a pedigree ( see VCFinjectPedigree ) or a pedigree must be defined.

## Lumpy-SV

 * 20180115: this tools recognize lumpy-sv genotypes


### see also

 *  VcfBurdenMAF
 *  VcfBurdenFilterExac

END_DOC
*/

@Program(name="vcfburdenfisherh",
	description="Fisher Case /Controls per Variant",
	keywords= {"vcf","burden","fisher"},
	modificationDate="20190418"
	)
public class VcfBurdenFisherH
	extends Launcher
	{	
	public static final double DEFAULT_MIN_FISHER_PVALUE = 0.05;
	private static final Logger LOG = Logger.build(VcfBurdenFisherH.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-fisher","--minFisherPValue"},description="if p-value fisher(case/control vs have alt/have not alt) lower than 'fisher' the FILTER Column is Filled")
	private double minFisherPValue = DEFAULT_MIN_FISHER_PVALUE ;
	@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Don't try to calculate things why variants already FILTERed (faster)")
	private boolean ignoreFiltered=false;
	@Parameter(names={"-p","--pedigree"},description="[20180115] Pedigree file. Default: use the pedigree data in the VCF header." + Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile=null;
	@Parameter(names={"-gtf","--gtf","--gtFiltered"},description="[20180115] Ignore FILTERed **Genotype**")
	private boolean ignore_filtered_genotype=false;
	@Parameter(names={"-lumpy-su-min","--lumpy-su-min"},description="[20180115] if variant identified as LUMPy-SV variant. This is the minimal number of 'SU' to consider the genotype as a variant.")
	private int lumpy_SU_threshold=1;
	@Parameter(names={"--attribute"},description="[20190418] Name of the attribue used as FILTER and INFO")
	private String burdenHFisherAttr = "BurdenHFisher";
	@Parameter(names={"--report"},description="[20190418] save report as bed file")
	private Path bedExportPath = null;

	
	

	
	private static class Count {
		int case_have_alt =0;
		int case_miss_alt = 0;
		int ctrl_have_alt = 0;
		int ctrl_miss_alt = 0;
		}
		
	public VcfBurdenFisherH() {
	}
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter w) {
		final Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor = 
				(header)->  new Pedigree.CaseControlExtractor().extract(header);
		
		final PrintWriter report;		
				
		if(this.bedExportPath==null) {
			report = new PrintWriter(new NullOuputStream());
			} 
		else
			{
			report = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(this.bedExportPath));
			}
		
		final Set<Pedigree.Person> individualSet;
		final VCFHeader header = r.getHeader();
		final boolean is_lumpy_vcf_header = LumpyConstants.isLumpyHeader(header);

		final VCFHeader h2= new VCFHeader(header);
		if(this.pedigreeFile == null)
			{
			individualSet = caseControlExtractor.apply(header);
			}
		else
			{
			try {
				individualSet = new Pedigree.CaseControlExtractor().extract(
						header,
						new Pedigree.Parser().parse(this.pedigreeFile)
						);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			}
		
		

		final VCFInfoHeaderLine fisherAlleleInfoHeader = new VCFInfoHeaderLine(
				this.burdenHFisherAttr ,
				VCFHeaderLineCount.A,
				VCFHeaderLineType.Float,
				"Fisher Exact Test Case/Control."
				);
		final VCFFilterHeaderLine fisherAlleleFilterHeader = new VCFFilterHeaderLine(
				fisherAlleleInfoHeader.getID(),
				"Fisher case:control vs miss|have ALT is lower than "+ this.minFisherPValue
				);
		
		final VCFInfoHeaderLine fisherDetailInfoHeader = new VCFInfoHeaderLine(
				this.burdenHFisherAttr + "Detail",
				VCFHeaderLineCount.A,VCFHeaderLineType.String,
				"Fisher Exact Test Case/Control"
				);
		h2.addMetaDataLine(fisherAlleleInfoHeader);
		h2.addMetaDataLine(fisherAlleleFilterHeader);
		h2.addMetaDataLine(fisherDetailInfoHeader);

		w.writeHeader(h2);
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(r.getHeader()).logger(LOG).build();
		while(r.hasNext())
			{
			final VariantContext ctx = progress.apply(r.next());
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.rmAttribute(fisherAlleleInfoHeader.getID());
			vcb.rmAttribute(fisherDetailInfoHeader.getID());
			
			final Set<String> oldFilters = new HashSet<>(ctx.getFilters());
			oldFilters.remove(fisherAlleleFilterHeader.getID());
			
			
			if(this.ignoreFiltered && ctx.isFiltered())
				{
				w.add(vcb.make());
				continue;
				}
			final boolean identified_as_lumpy= 
					is_lumpy_vcf_header && 
					LumpyConstants.isLumpyVariant(ctx)
					;
			boolean set_filter = true;
			boolean found_one_alt_to_compute = false;
			final List<String> infoData = new ArrayList<>(ctx.getAlleles().size());
			final List<Double> fisherValues = new ArrayList<>(ctx.getAlleles().size());
			
			for(final Allele observed_alt: ctx.getAlternateAlleles()) {
				if(observed_alt.isNoCall()) {
					infoData.add(
							String.join("|",
							"ALLELE",String.valueOf(observed_alt.getDisplayString()),
							"FISHER","-1.0"
							));
					fisherValues.add(-1.0);
					continue;
					}
				
				/* count for fisher allele */
				final Count count = new Count();
				
				/* loop over persons in this pop */
				for(final Pedigree.Person p: individualSet ) 	{
					/* get genotype for this individual */
					final Genotype genotype = ctx.getGenotype(p.getId());
					
					final boolean genotype_contains_allele;
					
					if(identified_as_lumpy && !genotype.isCalled())
						{
						if(this.ignore_filtered_genotype && genotype.isFiltered())
							{
							if(p.isAffected()) { count.case_miss_alt++; }
							else { count.ctrl_miss_alt++; }
							continue;
							}
						if(!genotype.hasExtendedAttribute("SU"))
							{
							throw new JvarkitException.FileFormatError(
									"Variant identified as lumpysv, but not attribute 'SU' defined in genotye "+genotype);
							}
						final int su_count = genotype.getAttributeAsInt("SU", 0);
						genotype_contains_allele = su_count>= this.lumpy_SU_threshold;
						}
					else
						{
						/* individual is not in vcf header */
						if(genotype==null || !genotype.isCalled() || (this.ignore_filtered_genotype && genotype.isFiltered())) {
							if(genotype==null) LOG.warn("Genotype is null for sample "+p.getId()+" not is pedigree!");
							//no information , we consider that sample was called AND HOM REF
							if(p.isAffected()) { count.case_miss_alt++; }
							else { count.ctrl_miss_alt++; }
							continue;
							}
						
						/* loop over alleles */
						genotype_contains_allele = genotype.getAlleles().stream().
								anyMatch(A->A.equals(observed_alt));
						
						}
					
					/* fisher */
					if(genotype_contains_allele) {
						if(p.isAffected()) { count.case_have_alt++; ;}
						else { count.ctrl_have_alt++; }
						}
					else {
						if(p.isAffected()) { count.case_miss_alt++; }
						else { count.ctrl_miss_alt++; }
						}
				}/* end of loop over persons */
				
				

				
				/* fisher test for alleles */
				final FisherExactTest fisherAlt = FisherExactTest.compute(
						count.case_have_alt, count.case_miss_alt,
						count.ctrl_have_alt, count.ctrl_miss_alt
						);
				
				fisherValues.add(fisherAlt.getAsDouble());
				infoData.add(
						String.join("|",
						"ALLELE",String.valueOf(observed_alt.getDisplayString()),
						"FISHER",String.valueOf(fisherAlt.getAsDouble()),
						"CASE_HAVE_ALT",String.valueOf(count.case_have_alt),
						"CASE_MISS_ALT",String.valueOf(count.case_miss_alt),
						"CTRL_HAVE_ALT",String.valueOf(count.ctrl_have_alt),
						"CTRL_MISS_ALT",String.valueOf(count.ctrl_miss_alt)
						));
				
				found_one_alt_to_compute = true;
				
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
				report.print(count.case_have_alt);
				report.print('\t');
				report.print(count.case_miss_alt);
				report.print('\t');
				report.print(count.ctrl_have_alt);
				report.print('\t');
				report.print(count.ctrl_miss_alt);
				report.print('\t');
				
				if( fisherAlt.getAsDouble() >= this.minFisherPValue ) {
					set_filter = false;
					report.print(fisherAlleleFilterHeader.getID());
					}
				else
					{
					report.print(".");
					}
				
				report.println();
				} //end of for each ALT allele

			vcb.attribute(fisherAlleleInfoHeader.getID(),fisherValues);
			vcb.attribute(fisherDetailInfoHeader.getID(),infoData );
			
			if( set_filter && found_one_alt_to_compute) {
				vcb.filter(fisherAlleleFilterHeader.getID());
				}
			w.add(vcb.make());
			report.println();
			}
		progress.close();
		w.close();
		report.flush();
		report.close();
		return 0;
		}

	
	@Override
	public int doWork(final List<String> args) {
		if(StringUtils.isBlank(this.burdenHFisherAttr)) {
			LOG.error("empty  burdenHFisherAttr");
			return -1;
			}
		return doVcfToVcf(args,this.outputFile);
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfBurdenFisherH().instanceMainWithExit(args);
		}
	}
