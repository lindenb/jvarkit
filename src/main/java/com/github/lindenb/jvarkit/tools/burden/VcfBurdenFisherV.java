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
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFBuffer;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;

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

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
		{
		@Parameter(names={"-if","--ignorefilter"},
			description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
		private boolean acceptFiltered = false;



		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final File tmpDir;
			private VCFBuffer tmpw = null;
			private final Map<Pedigree.Person,SuperVariant> indi2supervariant = new HashMap<>();
			private Count count= null;
			private VCFHeader header2 = null;
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				this.tmpDir = IOUtils.getDefaultTmpDir();
				}
			
			@Override
			public void writeHeader(final VCFHeader header) {
				this.indi2supervariant.clear();
				for(final Pedigree.Person  person: new Pedigree.CaseControlExtractor().extract(header)) {
					this.indi2supervariant.put(person, SuperVariant.SV0);
					}
				this.tmpw = new VCFBuffer(1000,tmpDir);
				this.tmpw.writeHeader(header);
				this.count = new Count();
				this.header2 = new VCFHeader(header);
				if(this.header2.getMetaDataLine(VCF_HEADER_FISHER_VALUE)!=null)
					{
					throw new JvarkitException.UserError(
							"VCF Header "+VCF_HEADER_FISHER_VALUE+" already specified in input");
					}
				}
			
			@Override
			public void add(final VariantContext ctx) {
				this.tmpw.add(ctx);
				
				if(ctx.isFiltered() && !CtxWriterFactory.this.acceptFiltered) return;
				final int n_alts = ctx.getAlternateAlleles().size();
				
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					return;
				}
				
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. Using getAltAlleleWithHighestAlleleCount.");
					}
				
				final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
				
				//loop over person in this pedigree
				for(final Pedigree.Person person : indi2supervariant.keySet() ) {
					if(this.indi2supervariant.get(person)==SuperVariant.AT_LEAST_ONE_VARIANT) continue;
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
							this.indi2supervariant.put(person,SuperVariant.AT_LEAST_ONE_VARIANT);
							break;
							}
						}//end of allele
					}//en dof for[person]			
				}
			
			@Override
			public void close() {
				VCFIterator in2  = null;
				try {
					for(final Pedigree.Person person : this.indi2supervariant.keySet() ) {
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
					
					this.header2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE,
							String.valueOf(fisher.getAsDouble())));
					this.header2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE+".count",
							String.join("|",
							"CASE_SV0="+count.count_case_sv0,
							"CASE_SV1="+count.count_case_sv1,
							"CTRL_SV0="+count.count_ctrl_sv0,
							"CTRL_SV1="+count.count_ctrl_sv1		
							)));
	
					in2 = this.tmpw.iterator();
					super.writeHeader(this.header2);
					while(in2.hasNext()) {
						super.add(in2.next());
						}
					}
				finally {
					CloserUtil.close(this.tmpw);
					if(this.tmpw!=null) this.tmpw.dispose();
					this.tmpw =null;
					CloserUtil.close(in2);
					this.indi2supervariant.clear();
					this.count=null;
					super.close();
					}
				}
			}
		
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		}
	
	
	public VcfBurdenFisherV()
		{
		}
	 
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter delegate)
		{
		final VariantContextWriter  out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
		out.writeHeader(in.getHeader());
		while(in.hasNext())
			{
			out.add(progess.watch(in.next()));
			}
		progess.finish();
		out.close();
		return 0;
		}
	
	@Override
	protected int doVcfToVcf(final List<String> args,final File outorNull) {
		return doVcfToVcfMultipleStream(oneFileOrNull(args), outorNull);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.component.initialize()!=0) {
				return -1;
				}
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	public static void main(String[] args)
		{
		new VcfBurdenFisherV().instanceMainWithExit(args);
		}
	}
