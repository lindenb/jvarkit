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

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/**

BEGIN_DOC


Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


### Output




#### INFO column


 *  BurdenMAFCas : MAF cases
 *  BurdenMAFControls : MAF controls





#### FILTER column


 *  BurdenMAFCas : MAF for cases  doesn't meet  user's requirements
 *  BurdenMAFControls : MAF for controls  doesn't meet  user's requirements
 *  BurdenMAFCaseOrControls : MAF for controls or cases  doesn't meet  user's requirements





### see also


 *  VcfBurdenFilterExac
 *  VcfBurdenFisherH


END_DOC
*/
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

Variant in that VCF should have one and **only one** ALT allele. Use [https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele](https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele) if needed.

### Output


#### INFO column

  * **BurdenMAFCas** : MAF cases
  * **BurdenMAFControls** : MAF controls

#### FILTER column

  * **BurdenMAFCas** : MAF for cases  doesn't meet  user's requirements
  * **BurdenMAFControls** : MAF for controls  doesn't meet  user's requirements
  * **BurdenMAFCaseOrControls** : MAF for controls or cases  doesn't meet  user's requirements

END_DOC

*/
@Program(name="vcfburdenmaf",
	description="Burden : MAF for Cases / Controls ",
	keywords={"vcf","burden","maf","case","control"}
	)
public class VcfBurdenMAF
	extends Launcher
	{
	private static final int CASE_POP=0;
	
	private static final Logger LOG = Logger.build(VcfBurdenMAF.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@Parameter(names={"-maxMAF","--maxMAF"},description="if MAF of cases OR MAF of control is greater than maxMAF, the the FILTER Column is Filled")
	private double maxMAF = 0.05 ;

	@Parameter(names={"-c","--homref"},description="Treat No Call './.' genotypes as HomRef")
	private boolean noCallAreHomRef = false;
	
	
	public VcfBurdenMAF()
		{
		}
	 
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final VCFHeader header=in.getHeader();
		final Pedigree pedigree = Pedigree.newParser().parse(header);
		if(pedigree.isEmpty())
			{
			LOG.error("No pedigree found in header "+inputName+". use VcfInjectPedigree to add it");
			return -1;
			}
		
		final Set<Pedigree.Person> caseSamples = pedigree.getAffected();
		final Set<Pedigree.Person> controlSamples = pedigree.getUnaffected();
		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));

			final VCFInfoHeaderLine mafCasInfoHeader = new VCFInfoHeaderLine(
					"BurdenMAFCas",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Cases"
					);
			final VCFInfoHeaderLine mafControlsInfoHeader = new VCFInfoHeaderLine(
					"BurdenMAFControls",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Controls"
					);
			final VCFFilterHeaderLine filterCasHeader = new VCFFilterHeaderLine(
					mafCasInfoHeader.getID(),"MAF of cases is greater than "+this.maxMAF
					);
			final VCFFilterHeaderLine filterControlsHeader = new VCFFilterHeaderLine(
					mafControlsInfoHeader.getID(),"MAF of controls is greater than "+this.maxMAF
					);
			final VCFFilterHeaderLine filterCaseOrControlsHeader = new VCFFilterHeaderLine(
					"BurdenMAFCaseOrControls","MAF of (cases OR controls) is greater than "+this.maxMAF
					);			
			
			h2.addMetaDataLine(mafCasInfoHeader);
			h2.addMetaDataLine(mafControlsInfoHeader);
			h2.addMetaDataLine(filterCasHeader);
			h2.addMetaDataLine(filterControlsHeader);
			h2.addMetaDataLine(filterCaseOrControlsHeader);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
			out.writeHeader(h2);
			while(in.hasNext())
				{
				final VariantContext ctx = progess.watch(in.next());
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				final List<Double> mafCasList = new ArrayList<>(); 
				final List<Double> mafCtrlList = new ArrayList<>(); 
				boolean set_max_maf_cas=true;
				boolean set_max_maf_control=true;
				boolean seen_data=false;
				
				for(final Allele observed_alt : ctx.getAlternateAlleles() )
					{
					/* loop over two populations : 0 = case, 1=controls */
					for(int pop=0;pop<2;++pop) {
						final MafCalculator mafCalculator = new MafCalculator(observed_alt, ctx.getContig());
						mafCalculator.setNoCallIsHomRef(this.noCallAreHomRef);
						
						
						/* loop over persons in this pop */
						for(final Pedigree.Person p:(pop==CASE_POP?caseSamples:controlSamples)) 
							{
							/* get genotype for this individual */
							final Genotype genotype = ctx.getGenotype(p.getId());
							mafCalculator.add(genotype, p.isMale());
							/* if(pop==CASE_POP && genotype.isCalled()) LOG.info("DEBUGMAF: "+p+" "+genotype); */
							}/* end of loop over persons */
						/* at least one genotype found */
						if(!mafCalculator.isEmpty())
							{
							seen_data=true;
							
							/* get MAF */
							final double maf = mafCalculator.getMaf();
							
							if(pop == CASE_POP) {
								/* add INFO attribute */
								mafCasList.add(maf);
								/* remove FILTER if needed */
								if(maf<=this.maxMAF)  set_max_maf_cas=false;
								}
							else
								{
								/* add INFO attribute */
								mafCtrlList.add(maf);
								/* remove FILTER if needed */
								if(maf<=this.maxMAF)  set_max_maf_control=false;
								}
							} 
						else
							{
							if(pop == CASE_POP) {
								mafCasList.add(-1.0);
								set_max_maf_cas=false;
							} else
							{
								mafCtrlList.add(-1.0);
								set_max_maf_control=false;
							}
							}
						}/* end of loop over pop */
					}/* end loop over alt allele */
				
				
				vcb.attribute(mafCasInfoHeader.getID(),mafCasList);
				vcb.attribute(mafControlsInfoHeader.getID(),mafCtrlList);
				
				if(seen_data) {
					if(set_max_maf_cas) vcb.filter(filterCasHeader.getID());
					if(set_max_maf_control) vcb.filter(filterControlsHeader.getID());
					if(set_max_maf_cas || set_max_maf_control) {
						vcb.filter(filterCaseOrControlsHeader.getID());
					}
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
		return doVcfToVcf(args,outputFile);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenMAF().instanceMainWithExit(args);
		}
	}
