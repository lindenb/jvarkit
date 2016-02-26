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
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.tools.vcfcmp.EqualRangeVcfIterator;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfBurdenFilter3
	extends AbstractVcfBurdenFilter3
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenFilter3.class);
	
	public VcfBurdenFilter3()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.exacFile==null || !super.exacFile.exists())
			{
			return wrapException("Undefined Exac file option -"+OPTION_EXACFILE);
			}
		return super.initializeKnime();
	 	}
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator vcfIterator,
			final VariantContextWriter out
			) throws IOException {
		VcfIterator exacIn =null;
		final VcfIterator in = VCFUtils.createAssertSortedVcfIterator(vcfIterator, VCFUtils.createTidPosComparator(vcfIterator.getHeader().getSequenceDictionary()));
		EqualRangeVcfIterator equalRange = null;
		final String exacPopulations[]= super.exacPopulationStr.split("[,]+");
		try {
			LOG.info("open "+super.exacFile);
			exacIn = VCFUtils.createVcfIteratorFromFile(super.exacFile);
			equalRange = new EqualRangeVcfIterator(exacIn, VCFUtils.createTidPosComparator(exacIn.getHeader().getSequenceDictionary()));
			final VCFHeader header=in.getHeader();
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine("BurdenF3",
					"Freq:"+this.maxFreq+"% Pop:"+super.exacPopulationStr);
			final VCFInfoHeaderLine freqExacInfoHeader = new VCFInfoHeaderLine(
					"FreqExac",1,VCFHeaderLineType.Float,"Freq in Exac AC/AN"
					);
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			h2.addMetaDataLine(freqExacInfoHeader);
			for(final String pop:exacPopulations)
				{
				final VCFInfoHeaderLine ac = exacIn.getHeader().getInfoHeaderLine("AC_"+pop);
				if(ac!=null) h2.addMetaDataLine(ac);
				final VCFInfoHeaderLine an = exacIn.getHeader().getInfoHeaderLine("AN_"+pop);
				if(an!=null) h2.addMetaDataLine(an);
				}
			h2.addMetaDataLine(filter);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				boolean set_filter = false;
				final VariantContext ctx = progess.watch(in.next());
				if(	ctx.getAlternateAlleles().size()!=1) {
					return wrapException("Expected only one allele per variant. Please use ManyAlleletoOne.");
					}
				final Allele alt = ctx.getAlternateAllele(0);
				
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

				Optional<Float> freqInExac = Optional.empty();
				for(final VariantContext exacVariant : equalRange.next(ctx))
					{
					if(!exacVariant.getReference().equals(ctx.getReference())) continue;
					int exacAltIndex=-1;
					for(final Allele  exacAlt : exacVariant.getAlternateAlleles()) {
						++exacAltIndex;
						if(!exacAlt.equals(alt)) continue;
						
						for(final String pop:exacPopulations)
							{
							final List<Object> aclist = exacVariant.getAttributeAsList("AC_"+pop);
							if(aclist.isEmpty() || aclist.size()<=exacAltIndex) continue;
							final List<Object> anlist = exacVariant.getAttributeAsList("AN_"+pop);
							if(anlist.isEmpty() || anlist.size()<=exacAltIndex) continue;

							
							
							
							final float ac = Float.parseFloat( aclist.get(exacAltIndex).toString());
							final float an = Float.parseFloat( anlist.get(exacAltIndex).toString());
							vcb.attribute("AC_"+pop,ac);
							vcb.attribute("AN_"+pop,an);
							
							if(an <= 0f) continue;
							final float freq = ac/an;
							if( !freqInExac.isPresent() || freq > freqInExac.get()) {
								freqInExac = Optional.of(freq);
								}
							}
						}
					}

				if(!freqInExac.isPresent() && super.ifNotInExacThenDiscard) {
					set_filter = true;
				}
				
				if(freqInExac.isPresent() && freqInExac.get() > (super.maxFreq)) {
					set_filter = true;
				}
				
				if(freqInExac.isPresent()) {
					vcb.attribute(freqExacInfoHeader.getID(), freqInExac.get());
				}
				if(  set_filter ) {
					vcb.filter(filter.getID());
					}
				out.add(vcb.make());
				}
			progess.finish();
			return RETURN_OK;
			} catch(final Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(equalRange);
				CloserUtil.close(exacIn);
				CloserUtil.close(in);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenFilter3().instanceMainWithExit(args);
		}
	}
