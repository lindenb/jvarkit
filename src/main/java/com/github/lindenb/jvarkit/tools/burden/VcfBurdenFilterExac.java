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

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.vcfcmp.EqualRangeVcfIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC

Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


### Output

#### INFO column


 *  FreqExac : Exac frequency.
 *  AC_* and AN_*: Transpose original population data from original Exac file


#### FILTER column

 *  BurdenExac : if FreqExac doesn't fit the criteria maxFreq


### see also


 *  VcfBurdenMAF


END_DOC
*/
@Program(name="vcfburdenexac",
		description="Burden filter 3 - Exac",
		keywords={"vcf","burden","exac"}
		)
public class VcfBurdenFilterExac
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenFilterExac.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-exac","--exac"},description="Path to Exac VCF file. At the time of writing, you'd better use a normalized version of Exac (see https://github.com/lindenb/jvarkit/wiki/VCFFixIndels )")
	private File exacFile = null;

	@Parameter(names={"-d","--discardNotInExac"},description="if variant was not found in Exac, discard the variant (set the FILTER). Default: don't set the FILTER.")
	private boolean ifNotInExacThenDiscard = false;

	@Parameter(names={"-maxFreq","--maxFreq"},description="set FILTER if max(exac frequency in any pop) is greater than this value)")
	private double maxFreq = 0.001 ;

	@Parameter(names={"-pop","--population"},description="comma separated populations in exac")
	private String exacPopulationStr = "AFR,AMR,EAS,FIN,NFE,SAS";

	@Parameter(names={"-tabix","--tabix"},description="use tabix index for Exac it is present. Might speed up things if the number of variant is low.")
	private boolean useTabixIndex = false;
	
	public VcfBurdenFilterExac()
		{
		}
	 
	@Override
	protected int doVcfToVcf(final String inputName,final VcfIterator vcfIterator,final VariantContextWriter out) {
		
		VcfIterator exacIn =null;
		TabixVcfFileReader tabix=null;
		final VcfIterator in = VCFUtils.createAssertSortedVcfIterator(vcfIterator, VCFUtils.createTidPosComparator(vcfIterator.getHeader().getSequenceDictionary()));
		EqualRangeVcfIterator equalRange = null;
		final String exacPopulations[]= this.exacPopulationStr.split("[,]+");
		final VCFHeader exacHeader;
		try {
			LOG.info("open "+this.exacFile);
			if(this.useTabixIndex && VCFUtils.isTabixVcfFile(this.exacFile))
				{
				tabix = new TabixVcfFileReader(this.exacFile.getPath());
				exacIn = null;
				equalRange = null;
				exacHeader = tabix.getHeader();
				}
			else
				{
				tabix=null;
				exacIn = VCFUtils.createVcfIteratorFromFile(this.exacFile);
				equalRange = new EqualRangeVcfIterator(exacIn, VCFUtils.createTidPosComparator(exacIn.getHeader().getSequenceDictionary()));
				exacHeader = exacIn.getHeader();
				}
			
			
			final VCFHeader header=in.getHeader();
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine("BurdenExac",
					"Freq:"+this.maxFreq+" Pop:"+this.exacPopulationStr);
			final VCFInfoHeaderLine freqExacInfoHeader = new VCFInfoHeaderLine(
					"FreqExac",1,VCFHeaderLineType.Float,"Freq in Exac AC/AN"
					);
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			h2.addMetaDataLine(freqExacInfoHeader);
			for(final String pop:exacPopulations)
				{
				final VCFInfoHeaderLine ac = exacHeader.getInfoHeaderLine("AC_"+pop);
				if(ac!=null) h2.addMetaDataLine(ac);
				final VCFInfoHeaderLine an = exacHeader.getInfoHeaderLine("AN_"+pop);
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
					LOG.error("Expected only one allele per variant. Please use ManyAlleletoOne.");
					return -1;
					}
				final Allele alt = ctx.getAlternateAllele(0);
				
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

				Optional<Float> freqInExac = Optional.empty();
				
				/* list of variant found in exac */
				final List<VariantContext> exacListOfVariants;
				
				if( equalRange!=null)  {
					exacListOfVariants =  equalRange.next(ctx);
				} else
				{
					exacListOfVariants = new ArrayList<>();
					final Iterator<VariantContext> vit=tabix.iterator(ctx.getContig(),ctx.getStart(), ctx.getEnd());
					while(vit.hasNext()) {
						final VariantContext exacv = vit.next();
						if(!exacv.getContig().equals(ctx.getContig())) continue;
						if(!exacv.getReference().equals(ctx.getReference())) continue;
						if(exacv.getStart()!=ctx.getStart()) continue;
						exacListOfVariants.add(exacv); 
					}
				}
				
				for(final VariantContext exacVariant :exacListOfVariants)
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

				if(!freqInExac.isPresent() && this.ifNotInExacThenDiscard) {
					set_filter = true;
				}
				
				if(freqInExac.isPresent() && freqInExac.get() > (this.maxFreq)) {
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
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(equalRange);
				CloserUtil.close(exacIn);
				CloserUtil.close(in);
				CloserUtil.close(tabix);
			}
		}
	
	 
	@Override
	public int doWork(List<String> args) {
		if(this.exacFile==null || !this.exacFile.exists())
			{
			LOG.error("Undefined Exac file option");
			return  -1;
			}
		return doVcfToVcf(args, outputFile);
		}
	
	
	public static void main(String[] args)
		{
		new VcfBurdenFilterExac().instanceMainWithExit(args);
		}
	}
