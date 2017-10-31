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
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

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
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
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
import com.beust.jcommander.ParametersDelegate;
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

	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	@XmlType(name="vcfburdenmaf")
	@XmlRootElement(name="vcfburdenmaf")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
		{
		@XmlElement(name="max-maf")
		@Parameter(names={"-maxMAF","--maxMAF"},description="if MAF of cases OR MAF of control is greater than maxMAF, the the FILTER Column is Filled")
		private double maxMAF = 0.05 ;
		
		@XmlElement(name="nocall-is-homref")
		@Parameter(names={"-c","--homref"},description="Treat No Call './.' genotypes as HomRef")
		private boolean noCallAreHomRef = false;
		
		@XmlElement(name="ignore-filtered")
		@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Don't try to calculate things why variants already FILTERed (faster)")
		private boolean ignoreFiltered=false;
		
		@XmlTransient
		private Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor = 
			(header)->  new Pedigree.CaseControlExtractor().extract(header);
		
		/** supplier for a collection of case/control individual.
		 * default will use the pedigree. The idea
		 * is to get a functional interface to give a chance to provide a
		 * list of sample when using vcfburdenfisherh as a java bean 
		 */
		public void setCaseControlExtractor(final Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor) {
			this.caseControlExtractor = caseControlExtractor;
		}
		
		public void setMaxMAF(double maxMAF) {
			this.maxMAF = maxMAF;
		}
		
		public void setNoCallAreHomRef(boolean noCallAreHomRef) {
			this.noCallAreHomRef = noCallAreHomRef;
		}
		
		public void setIgnoreFiltered(boolean ignoreFiltered) {
			this.ignoreFiltered = ignoreFiltered;
		}

		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final boolean ignoreFiltered = CtxWriterFactory.this.ignoreFiltered;
			private final Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor = CtxWriterFactory.this.caseControlExtractor;
			private Set<Pedigree.Person> caseSamples=null;
			private Set<Pedigree.Person> controlSamples=null;
			private final VCFInfoHeaderLine mafCasInfoHeader = new VCFInfoHeaderLine(
					"BurdenMAFCas",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Cases"
					);
			private final VCFInfoHeaderLine mafControlsInfoHeader = new VCFInfoHeaderLine(
					"BurdenMAFControls",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Controls"
					);
			private final VCFFilterHeaderLine filterCasHeader = new VCFFilterHeaderLine(
					mafCasInfoHeader.getID(),"MAF of cases is greater than "+CtxWriterFactory.this.maxMAF
					);
			private final VCFFilterHeaderLine filterControlsHeader = new VCFFilterHeaderLine(
					mafControlsInfoHeader.getID(),"MAF of controls is greater than "+CtxWriterFactory.this.maxMAF
					);
			private final VCFFilterHeaderLine filterCaseOrControlsHeader = new VCFFilterHeaderLine(
					"BurdenMAFCaseOrControls","MAF of (cases OR controls) is greater than "+CtxWriterFactory.this.maxMAF
					);			

			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			@Override
			public void writeHeader(final VCFHeader header) {				
				final Set<Pedigree.Person> persons = this.caseControlExtractor.apply(header);
				
				this.caseSamples = persons.stream().
						filter(I->I.isAffected()).collect(Collectors.toSet());
				this.controlSamples = persons.stream().
						filter(I->I.isUnaffected()).collect(Collectors.toSet());

				final VCFHeader h2= new VCFHeader(header);
				h2.addMetaDataLine(this.mafCasInfoHeader);
				h2.addMetaDataLine(this.mafControlsInfoHeader);
				h2.addMetaDataLine(this.filterCasHeader);
				h2.addMetaDataLine(this.filterControlsHeader);
				h2.addMetaDataLine(this.filterCaseOrControlsHeader);
				super.writeHeader(h2);
				}
			
			@Override
			public void add(final VariantContext ctx) {
				if(this.ignoreFiltered && ctx.isFiltered())
					{
					super.add(ctx);
					return;
					}
				
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
						mafCalculator.setNoCallIsHomRef(CtxWriterFactory.this.noCallAreHomRef);
						
						/* loop over persons in this pop */
						for(final Pedigree.Person p:(pop==CASE_POP?this.caseSamples:this.controlSamples)) 
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
								if(maf<=CtxWriterFactory.this.maxMAF)  set_max_maf_cas=false;
								}
							else
								{
								/* add INFO attribute */
								mafCtrlList.add(maf);
								/* remove FILTER if needed */
								if(maf<=CtxWriterFactory.this.maxMAF)  set_max_maf_control=false;
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
				
				
				vcb.attribute(this.mafCasInfoHeader.getID(),mafCasList);
				vcb.attribute(this.mafControlsInfoHeader.getID(),mafCtrlList);
				
				if(seen_data) {
					if(set_max_maf_cas) vcb.filter(this.filterCasHeader.getID());
					if(set_max_maf_control) vcb.filter(this.filterControlsHeader.getID());
					if(set_max_maf_cas || set_max_maf_control) {
						vcb.filter(this.filterCaseOrControlsHeader.getID());
						}
					}

				super.add(vcb.make());
				}
			@Override
			public void close() {
				super.close();
				this.caseSamples=null;
				this.controlSamples=null;
				}
			}
		@Override
		public VariantContextWriter open(final VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		}
	
	
	public VcfBurdenMAF()
		{
		}
	 	
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VcfIterator in,
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
		new VcfBurdenMAF().instanceMainWithExit(args);
		}
	}
