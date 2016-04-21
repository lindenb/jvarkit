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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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

/**
 *     
 * VcfBurdenFisherH
 * @author lindenb
 *
 */
public class VcfBurdenFisherH
	extends AbstractVcfBurdenFisherH
	{	
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenFisherH.class);
	
	private static class Count {
		int case_have_alt =0;
		int case_miss_alt = 0;
		int ctrl_have_alt = 0;
		int ctrl_miss_alt = 0;
	}
	
	public VcfBurdenFisherH()
		{
		}
	 
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		final VCFHeader header=in.getHeader();
		final Set<Pedigree.Person> individuals = super.getCasesControlsInPedigree(header);
		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			
			final VCFInfoHeaderLine fisherAlleleInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisher",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Fisher Exact Test Case/Control."
					);
			final VCFFilterHeaderLine fisherAlleleFilterHeader = new VCFFilterHeaderLine(
					fisherAlleleInfoHeader.getID(),"Fisher case:control vs miss|have ALT is lower than "+super.minFisherPValue
					);
			
			final VCFInfoHeaderLine fisherDetailInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisherDetail",VCFHeaderLineCount.A,VCFHeaderLineType.String,"Fisher Exact Test Case/Control"
					);

			
			h2.addMetaDataLine(fisherAlleleInfoHeader);
			h2.addMetaDataLine(fisherAlleleFilterHeader);
			h2.addMetaDataLine(fisherDetailInfoHeader);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
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
					if( fisherAlt.getAsDouble() >= super.minFisherPValue ) {
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
				return wrapException(err);
			} finally {
				CloserUtil.close(in);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenFisherH().instanceMainWithExit(args);
		}
	}
