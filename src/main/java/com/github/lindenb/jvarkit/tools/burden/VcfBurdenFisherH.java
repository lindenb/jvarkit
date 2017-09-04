/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
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

import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC


Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.
VCF header must contain a pedigree ( see VCFinjectPedigree ).




### Output




#### INFO column


 *  BurdenFisher : Fisher test





#### FILTER column


 *  BurdenFisher :Fisher test doesn't meet  user's requirements





### see also


 *  VcfBurdenMAF
 *  VcfBurdenFilterExac






END_DOC
*/

@Program(name="vcfburdenfisherh",description="Fisher Case /Controls per Variant")
public class VcfBurdenFisherH
	extends Launcher
	{	
	private static final Logger LOG = Logger.build(VcfBurdenFisherH.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@Parameter(names={"-fisher","--minFisherPValue"},description="if p-value fisher(case/control vs have alt/have not alt) lower than 'fisher' the FILTER Column is Filled")
	private double minFisherPValue = 0.05 ;
	
	private static class Count {
		int case_have_alt =0;
		int case_miss_alt = 0;
		int ctrl_have_alt = 0;
		int ctrl_miss_alt = 0;
	}
	
	public VcfBurdenFisherH()
		{
		}
	 
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final VCFHeader header=in.getHeader();
		final Set<Pedigree.Person> individuals = super.getCasesControlsInPedigree(header);
		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			
			final VCFInfoHeaderLine fisherAlleleInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisher",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Fisher Exact Test Case/Control."
					);
			final VCFFilterHeaderLine fisherAlleleFilterHeader = new VCFFilterHeaderLine(
					fisherAlleleInfoHeader.getID(),"Fisher case:control vs miss|have ALT is lower than "+this.minFisherPValue
					);
			
			final VCFInfoHeaderLine fisherDetailInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisherDetail",VCFHeaderLineCount.A,VCFHeaderLineType.String,"Fisher Exact Test Case/Control"
					);

			
			h2.addMetaDataLine(fisherAlleleInfoHeader);
			h2.addMetaDataLine(fisherAlleleFilterHeader);
			h2.addMetaDataLine(fisherDetailInfoHeader);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
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
					for(final Pedigree.Person p:individuals ) 	{
						/* get genotype for this individual */
						final Genotype genotype = ctx.getGenotype(p.getId());
						/* individual is not in vcf header */
						if(genotype==null || !genotype.isCalled() ) {
							//no information , we consider that sample was called AND HOM REF
							if(p.isAffected()) { count.case_miss_alt++; }
							else { count.ctrl_miss_alt++; }
							continue;
						}
						
						/* loop over alleles */
						final boolean genotype_contains_allele=genotype.getAlleles().stream().anyMatch(A->A.equals(observed_alt));
						
						
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
					if( fisherAlt.getAsDouble() >= this.minFisherPValue ) {
						set_filter = false;
						}
					} //end of for each ALT allele
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

				vcb.attribute(fisherAlleleInfoHeader.getID(),fisherValues);
				vcb.attribute(fisherDetailInfoHeader.getID(),infoData );
				
				if( set_filter && found_one_alt_to_compute) {
					vcb.filter(fisherAlleleFilterHeader.getID());
					}
				
				out.add(vcb.make());
				}
			progess.finish();
			LOG.info("done");
			return RETURN_OK;
			} catch(Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(in);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args, outputFile);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenFisherH().instanceMainWithExit(args);
		}
	}
