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
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.vcfcmp.EqualRangeVcfIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC

20170626: this tool now supports multiple ALT in the user VCF, however it's not been tested for choosing when to set the FILTER or the min value

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

	@Parameter(names={"-d","--discardNotInExac"},description="if variant was not found in Exac, set the FILTER. Default: don't set the FILTER.")
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
		VCFFileReader tabix=null;
		final VcfIterator in = VCFUtils.createAssertSortedVcfIterator(vcfIterator, VCFUtils.createTidPosComparator(vcfIterator.getHeader().getSequenceDictionary()));
		EqualRangeVcfIterator equalRange = null;
		final String exacPopulations[]= this.exacPopulationStr.split("[,]+");
		final VCFHeader exacHeader;
		try {
			LOG.info("open "+this.exacFile);
			if(this.useTabixIndex && VCFUtils.isTabixVcfFile(this.exacFile))
				{
				tabix = new VCFFileReader(this.exacFile,true);
				exacIn = null;
				equalRange = null;
				exacHeader = tabix.getFileHeader();
				}
			else
				{
				tabix=null;
				exacIn = VCFUtils.createVcfIteratorFromFile(this.exacFile);
				equalRange = new EqualRangeVcfIterator(exacIn, VCFUtils.createTidPosComparator(exacIn.getHeader().getSequenceDictionary()));
				exacHeader = exacIn.getHeader();
				}
			
			
			final VCFHeader header=in.getHeader();
			final VCFHeader h2=addMetaData(new VCFHeader(header));
			

			
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine("BurdenExac",
					"Freq:"+this.maxFreq+" Pop:"+this.exacPopulationStr);
			h2.addMetaDataLine(filter);
			
			final VCFInfoHeaderLine freqExacInfoHeader = new VCFInfoHeaderLine(
					"FreqExac",
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Freq in Exac AC/AN"
					);
			h2.addMetaDataLine(freqExacInfoHeader);
			
			
			for(final String pop:exacPopulations)
				{
				if(pop.isEmpty()) continue;
				final VCFInfoHeaderLine ac = exacHeader.getInfoHeaderLine("AC_"+pop);
				if(ac==null) 
					{
					LOG.error("Cannot find AC_"+pop+" in "+this.exacFile+" header");
					return -1;
					}
				if(ac.getCountType()!=VCFHeaderLineCount.A)
					{
					LOG.error("expected VCFHeaderLineCount.A for "+ac );
					return -1;
					}
				h2.addMetaDataLine(ac);
				final VCFInfoHeaderLine an = exacHeader.getInfoHeaderLine("AN_"+pop);
				if(an==null)
					{
					LOG.error("Cannot find AN_"+pop+" in "+this.exacFile+" header");
					return -1;
					}
				if(an.getCount()!=1)
					{
					LOG.error("expected getCount==1 for "+an );
					return -1;
					}
				h2.addMetaDataLine(an);
				}
			
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
			out.writeHeader(h2);
			while(in.hasNext())
				{
				boolean set_filter = false;
				final VariantContext ctx = progess.watch(in.next());
				final List<Allele> ctxAlternateAlleles = ctx.getAlternateAlleles();
				if(	ctxAlternateAlleles.size()!=1) {
					LOG.error("Warning found more than one ALT allele per variant ("+ctx.getContig()+":"+ctx.getStart()+"/"+ctx.getReference()+"). Please use ManyAlleletoOne if you want to split those alleles");
					return -1;
					}
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

				
				/* list of variant found in exac */
				final List<VariantContext> exacListOfVariants;
				
				if( equalRange!=null)  {
					exacListOfVariants =  equalRange.next(ctx);
				} else
				{
					exacListOfVariants = new ArrayList<>();
					final CloseableIterator<VariantContext> vit=tabix.query(ctx.getContig(),ctx.getStart(), ctx.getEnd());
					while(vit.hasNext()) {
						final VariantContext exacv = vit.next();
						if(!exacv.getContig().equals(ctx.getContig())) continue;
						if(!exacv.getReference().equals(ctx.getReference())) continue;
						if(exacv.getStart()!=ctx.getStart()) continue;
						exacListOfVariants.add(exacv); 
						}
					vit.close();
				}
				
				if(exacListOfVariants.size()>1)
					{
					LOG.warn("There is more than one variant for ("+ctx.getContig()+":"+ctx.getStart()+"/"+ctx.getReference()+") in "+exacFile);
					}
				
				// AF for each allele
				final List<Float> freqInExac = new ArrayList<>(ctxAlternateAlleles.size());
				while(freqInExac.size()< ctxAlternateAlleles.size())
					{
					freqInExac.add(null);
					}
				
				/* loop over populations */
				for(final String pop:exacPopulations)
					{
					/* will be the AC count for this population */
					final List<Integer> acAttribute= new ArrayList<>(ctxAlternateAlleles.size());
					/* will be the AN count for this population */
					Integer anAttribute= null;
					/* loop over user ALT alleles */
					for(int ctxAltIdx=0;ctxAltIdx <ctxAlternateAlleles.size();++ctxAltIdx )
						{
						final Allele ctxAlt = ctxAlternateAlleles.get(ctxAltIdx);
						Integer exacAC=null;
						
						/* loop over exac variants */
						for(final VariantContext exacVariant :exacListOfVariants)
							{
							if(!exacVariant.getReference().equals(ctx.getReference())) continue;
							final List<Object> aclist = exacVariant.getAttributeAsList("AC_"+pop);
							if(aclist.isEmpty()) continue;
							final List<Object> anlist = exacVariant.getAttributeAsList("AN_"+pop);
							if(anlist.isEmpty() || anlist.size()!=1) continue;
							final int exacan = Integer.parseInt( anlist.get(0).toString());
							
							final List<Allele> exacAlternateAlleles = exacVariant.getAlternateAlleles();
							
							/* loop over exac ALLeles */
							for(int exacAltIndex=0;exacAltIndex < exacAlternateAlleles.size();++exacAltIndex ) {
								final Allele exacAlt = exacAlternateAlleles.get(exacAltIndex);
								if(!exacAlt.equals(ctxAlt)) continue;
								if(aclist.isEmpty() || aclist.size()<=exacAltIndex) continue;

								final int ac = Integer.parseInt( aclist.get(exacAltIndex).toString());
								
								exacAC = ac;
								anAttribute = exacan;
								
								
								if(exacan <= 0) continue;
								final float freq = ac/exacan;
								if( freqInExac.get(ctxAltIdx)==null || freq > freqInExac.get(ctxAltIdx)) {
									freqInExac.set(ctxAltIdx, freq);
									}
								}
							}
						acAttribute.add(exacAC);
						}
					if(anAttribute!=null && acAttribute.stream().filter(N->N!=null).findAny().isPresent())
						{
						vcb.attribute("AC_"+pop,acAttribute);
						vcb.attribute("AN_"+pop,anAttribute);
						}
					}

				if(!freqInExac.stream().filter(F->F!=null).findAny().isPresent() && this.ifNotInExacThenDiscard) {
					set_filter = true;
				}
				
				if( freqInExac.stream().filter(F->F!=null).findAny().isPresent() && 
					freqInExac.stream().filter(F->F!=null).mapToDouble(F->F.doubleValue()).min().getAsDouble() > (this.maxFreq)) {
					set_filter = true;
				}
				
				if(freqInExac.stream().filter(F->F!=null).findAny().isPresent()) {
					vcb.attribute(freqExacInfoHeader.getID(), freqInExac);
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
