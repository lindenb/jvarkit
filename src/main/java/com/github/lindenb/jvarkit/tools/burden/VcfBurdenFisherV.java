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
import java.util.Collection;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 *     
 * VcfBurdenFisherV
 * @author lindenb
 *
 */
public class VcfBurdenFisherV
	extends AbstractVcfBurdenFisherV
	{
	public static final String VCF_HEADER_FISHER_VALUE="VCFBurdenFisherV";

	private enum SuperVariant
		{
		SV0,AT_LEAST_ONE_VARIANT
		}
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenFisherV.class);
	
	private static class Count {
		int count_case_sv0 =0;
		int count_ctrl_sv0 = 0;
		int count_case_sv1 = 0;
		int count_ctrl_sv1 = 0;
	}
	
	public VcfBurdenFisherV()
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
		final Set<Pedigree.Person> individuals =super.getCasesControlsInPedigree(header);
		
		
		
		VCFBuffer tmpw = null;
		VcfIterator in2=null;
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			tmpw = new VCFBuffer(1000,this.getTmpdir());
			tmpw.writeHeader(header);
			final Count count= new Count();
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			while(in.hasNext()) {
				final VariantContext ctx = progess.watch(in.next());
				tmpw.add(ctx);
				
				if(ctx.isFiltered() && !super.acceptFiltered) continue;
				int n_alts = ctx.getAlternateAlleles().size();
				
				
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					continue;
				}
				
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. Using getAltAlleleWithHighestAlleleCount.");
				}
				
				SuperVariant superVariant = SuperVariant.SV0;
				final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
				//loop over person in this pedigree
				for(final Pedigree.Person person : individuals) {
					final Genotype g = ctx.getGenotype(person.getId());	
					if(g==null) continue;//not in vcf header
					if(g.isFiltered()) continue;//ignore this genotype
					for(final Allele alt : g.getAlleles()) {
						if(observed_alt.equals(alt)) {
							superVariant = SuperVariant.AT_LEAST_ONE_VARIANT;
							break;
							}
						}//end of allele
					
					if(superVariant==SuperVariant.SV0 ) {
						if(person.isAffected()) count.count_case_sv0++;
						else count.count_ctrl_sv0++;
					} else // AT_LEAST_ONE_VARIANT 
						{
						if(person.isAffected()) count.count_case_sv1++;
						else count.count_ctrl_sv1++;
						}
					}//end of person
				}
			tmpw.close();
			progess.finish();
			
			
			final FisherExactTest fisher = FisherExactTest.compute(
					count.count_case_sv0, count.count_case_sv1,
					count.count_ctrl_sv0, count.count_ctrl_sv1
					);
			LOG.info("Fisher "+fisher.getAsDouble());
			h2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE,
					String.valueOf(fisher.getAsDouble())));
			h2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE+".count",
					String.join("|",
					"CASE_SV0="+count.count_case_sv0,
					"CASE_SV1="+count.count_case_sv1,
					"CTRL_SV0="+count.count_ctrl_sv0,
					"CTRL_SV1="+count.count_ctrl_sv1		
					)));

			in2 = tmpw.iterator();
			final SAMSequenceDictionaryProgress progess2=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in2.hasNext() &&  !out.checkError()) {
				final VariantContext ctx = progess2.watch(in2.next());
				out.add(ctx);
			}
			progess2.finish();
			
			LOG.info("done");
			return RETURN_OK;
			} catch(Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(tmpw);
				if(tmpw!=null) tmpw.dispose();
				CloserUtil.close(in);
				CloserUtil.close(in2);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args)
		{
		new VcfBurdenFisherV().instanceMainWithExit(args);
		}
	}
