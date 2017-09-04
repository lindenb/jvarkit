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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.beust.jcommander.Parameter;

/**

BEGIN_DOC

Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.

### Output


#### INFO column

 *  BurdenF1Fisher : Fisher test


#### FILTER column

 *  BurdenF1Fisher :Fisher test doesn't meet  user's requirements


### see also


 *  VcfBurdenFilter3


END_DOC
*/

@Program(
		name="vcfburdenfisherv",
		description="Fisher Case / Controls per Variant (Vertical)",
		keywords={"vcf","burden","fisher"}
		)
public class VcfBurdenFisherV
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenFisherV.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-if","--ignorefilter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;

	public static final String VCF_HEADER_FISHER_VALUE="VCFBurdenFisherV";

	private enum SuperVariant
		{
		SV0,AT_LEAST_ONE_VARIANT
		}
	
	private static class Count {
		int count_case_sv0 =0;
		int count_ctrl_sv0 = 0;
		int count_case_sv1 = 0;
		int count_ctrl_sv1 = 0;
	}
	
	public VcfBurdenFisherV()
		{
		}
	 
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final File tmpDir = new File(System.getProperty("java.io.tmpdir"));
		final VCFHeader header=in.getHeader();
		
		final Map<Pedigree.Person,SuperVariant> indi2supervariant = new HashMap<>();
		for(final Pedigree.Person  person: super.getCasesControlsInPedigree(header)) {
			indi2supervariant.put(person, SuperVariant.SV0);
		}
		
		VCFBuffer tmpw = null;
		VcfIterator in2=null;
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			tmpw = new VCFBuffer(1000,tmpDir);
			tmpw.writeHeader(header);
			final Count count= new Count();
		
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(in.hasNext()) {
				final VariantContext ctx = progess.watch(in.next());
				tmpw.add(ctx);
				
				if(ctx.isFiltered() && !this.acceptFiltered) continue;
				int n_alts = ctx.getAlternateAlleles().size();
				
				
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					continue;
				}
				
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. Using getAltAlleleWithHighestAlleleCount.");
				}
				
				final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
				
				//loop over person in this pedigree
				for(final Pedigree.Person person : indi2supervariant.keySet() ) {
					if(indi2supervariant.get(person)==SuperVariant.AT_LEAST_ONE_VARIANT) continue;
					final Genotype g = ctx.getGenotype(person.getId());	
					if(g==null) {
						continue;//not in vcf header
					}
					if(g.isFiltered()) {
						LOG.warn("ignoring filtered genotype");
						continue;//not filter.
					}
					for(final Allele alt : g.getAlleles()) {
						if(observed_alt.equals(alt)) {
							indi2supervariant.put(person,SuperVariant.AT_LEAST_ONE_VARIANT);
							break;
							}
						}//end of allele
					}//en dof for[person]
				} //end variant iteration
			tmpw.close();
			progess.finish();
			
			for(final Pedigree.Person person : indi2supervariant.keySet() ) {
					final SuperVariant superVariant = indi2supervariant.get(person);
					if(superVariant==SuperVariant.SV0 ) {
						if(person.isAffected()) count.count_case_sv0++;
						else count.count_ctrl_sv0++;
					} else // AT_LEAST_ONE_VARIANT 
						{
						if(person.isAffected()) count.count_case_sv1++;
						else count.count_ctrl_sv1++;
						}
				}//end of person
			
			
			
			
			final FisherExactTest fisher = FisherExactTest.compute(
					count.count_case_sv0, count.count_case_sv1,
					count.count_ctrl_sv0, count.count_ctrl_sv1
					);
			LOG.info("Fisher "+fisher.getAsDouble());
			if(h2.getMetaDataLine(VCF_HEADER_FISHER_VALUE)!=null)
				{
				LOG.error("VCF Header "+VCF_HEADER_FISHER_VALUE+" already specified in input");
				return -1;
				}
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
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(tmpw);
				if(tmpw!=null) tmpw.dispose();
				CloserUtil.close(in);
				CloserUtil.close(in2);
			}
		}
	
	@Override
	protected int doVcfToVcf(List<String> args, File outorNull) {
		return doVcfToVcfMultipleStream(oneFileOrNull(args), outorNull);
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args,this.outputFile);
		}
	
	public static void main(String[] args)
		{
		new VcfBurdenFisherV().instanceMainWithExit(args);
		}
	}
