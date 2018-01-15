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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import org.eclipse.jetty.io.RuntimeIOException;

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

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC


Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.
VCF header must contain a pedigree ( see VCFinjectPedigree ) or a pedigree must be defined.

## Lumpy-SV

 * 20180115: this tools recognize lumpy-sv genotypes


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

@Program(name="vcfburdenfisherh",
	description="Fisher Case /Controls per Variant",
	keywords= {"vcf","burden","fisher"})
public class VcfBurdenFisherH
	extends Launcher
	{	
	public static final double DEFAULT_MIN_FISHER_PVALUE = 0.05;
	private static final Logger LOG = Logger.build(VcfBurdenFisherH.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	@XmlType(name="vcfburdenfisherh")
	@XmlRootElement(name="vcfburdenfisherh")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		@XmlElement(name="min-fisher")
		@Parameter(names={"-fisher","--minFisherPValue"},description="if p-value fisher(case/control vs have alt/have not alt) lower than 'fisher' the FILTER Column is Filled")
		private double minFisherPValue = DEFAULT_MIN_FISHER_PVALUE ;
		@XmlElement(name="ignore-filtered")
		@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Don't try to calculate things why variants already FILTERed (faster)")
		private boolean ignoreFiltered=false;
		@XmlElement(name="pedigree")
		@Parameter(names={"-p","--pedigree"},description="[20180115] Pedigree file. Default: use the pedigree data in the VCF header." + Pedigree.OPT_DESCRIPTION)
		private File pedigreeFile=null;
		@XmlElement(name="ignore-filtered-genotype")
		@Parameter(names={"-gtf","--gtf","--gtFiltered"},description="[20180115] Ignore FILTERed **Genotype**")
		private boolean ignore_filtered_genotype=false;
		@XmlElement(name="lumpy-su-min")
		@Parameter(names={"-lumpy-su-min","--lumpy-su-min"},description="[20180115] if variant identified as LUMPy-SV variant. This is the minimal number of 'SU' to consider the genotype as a variant.")
		private int lumpy_SU_threshold=1;

		
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
		public void setIgnoreFiltered(boolean ignoreFiltered) {
			this.ignoreFiltered = ignoreFiltered;
		}
		public void setMinFisherPValue(double minFisherPValue) {
			this.minFisherPValue = minFisherPValue;
		}
		
		
		private static class Count {
			int case_have_alt =0;
			int case_miss_alt = 0;
			int ctrl_have_alt = 0;
			int ctrl_miss_alt = 0;
		}
		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final VCFInfoHeaderLine fisherAlleleInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisher",VCFHeaderLineCount.A,VCFHeaderLineType.Float,
					"Fisher Exact Test Case/Control."
					);
			private final VCFFilterHeaderLine fisherAlleleFilterHeader = new VCFFilterHeaderLine(
					fisherAlleleInfoHeader.getID(),
					"Fisher case:control vs miss|have ALT is lower than "+CtxWriterFactory.this.minFisherPValue
					);
			
			private final VCFInfoHeaderLine fisherDetailInfoHeader = new VCFInfoHeaderLine(
					"BurdenHFisherDetail",
					VCFHeaderLineCount.A,VCFHeaderLineType.String,
					"Fisher Exact Test Case/Control"
					);
			private Set<Pedigree.Person> individualSet= null;
			private final boolean ignoreFiltered = CtxWriterFactory.this.ignoreFiltered;
			private final Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor = CtxWriterFactory.this.caseControlExtractor;
			private final boolean ignore_filtered_genotype = CtxWriterFactory.this.ignore_filtered_genotype;
			private final int lumpy_SU_threshold =  CtxWriterFactory.this.lumpy_SU_threshold;
			
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
			}
			
			@Override
			public void writeHeader(final VCFHeader header) {
				final VCFHeader h2= new VCFHeader(header);
				h2.addMetaDataLine(this.fisherAlleleInfoHeader);
				h2.addMetaDataLine(this.fisherAlleleFilterHeader);
				h2.addMetaDataLine(this.fisherDetailInfoHeader);
				if( CtxWriterFactory.this.pedigreeFile == null)
					{
					this.individualSet = this.caseControlExtractor.apply(header);
					}
				else
					{
					try {
						this.individualSet = new Pedigree.CaseControlExtractor().extract(
								header,
								new Pedigree.Parser().parse( CtxWriterFactory.this.pedigreeFile)
								);
						}
					catch(final IOException err)
						{
						throw new RuntimeIOException(err);
						}
					}
				super.writeHeader(h2);
				}
			
			@Override
			public void add(final VariantContext ctx) {
				if(this.ignoreFiltered && ctx.isFiltered())
					{
					super.add(ctx);
					return;
					}
				final boolean identified_as_lumpy= 
						ctx.getStructuralVariantType()!=null &&
						ctx.getAlternateAlleles().size()==1 &&
						ctx.getAlternateAllele(0).isSymbolic() &&
						ctx.hasAttribute("SU")
						;
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
					for(final Pedigree.Person p: this.individualSet ) 	{
						/* get genotype for this individual */
						final Genotype genotype = ctx.getGenotype(p.getId());
						
						final boolean genotype_contains_allele;
						
						if(identified_as_lumpy && !genotype.isCalled())
							{
							if(this.ignore_filtered_genotype && genotype.isFiltered())
								{
								if(p.isAffected()) { count.case_miss_alt++; }
								else { count.ctrl_miss_alt++; }
								continue;
								}
							if(!genotype.hasExtendedAttribute("SU"))
								{
								throw new JvarkitException.FileFormatError(
										"Variant identified as lumpysv, but not attribute 'SU' defined in genotye "+genotype);
								}
							final int su_count = genotype.getAttributeAsInt("SU", 0);
							genotype_contains_allele = su_count>= this.lumpy_SU_threshold;
							}
						else
							{
							/* individual is not in vcf header */
							if(genotype==null || !genotype.isCalled() || (this.ignore_filtered_genotype && genotype.isFiltered())) {
								if(genotype==null) LOG.warn("Genotype is null for sample "+p.getId()+" not is pedigree!");
								//no information , we consider that sample was called AND HOM REF
								if(p.isAffected()) { count.case_miss_alt++; }
								else { count.ctrl_miss_alt++; }
								continue;
								}
							
							/* loop over alleles */
							genotype_contains_allele = genotype.getAlleles().stream().
									anyMatch(A->A.equals(observed_alt));
							
							}
						
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
					if( fisherAlt.getAsDouble() >= CtxWriterFactory.this.minFisherPValue ) {
						set_filter = false;
						}
					} //end of for each ALT allele
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);

				vcb.attribute(this.fisherAlleleInfoHeader.getID(),fisherValues);
				vcb.attribute(this.fisherDetailInfoHeader.getID(),infoData );
				
				if( set_filter && found_one_alt_to_compute) {
					vcb.filter(this.fisherAlleleFilterHeader.getID());
					}
				
				super.add(vcb.make());
				}
			@Override
			public void close() {
				super.close();
				this.individualSet=null;
				}
			}
		
		
		@Override
		public VariantContextWriter open(VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		}
	
	public VcfBurdenFisherH() {
	}
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VcfIterator r, final VariantContextWriter delegate) {
		final VariantContextWriter w = this.component.open(delegate);
		w.writeHeader(r.getHeader());
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		while(r.hasNext())
			{
			w.add(progress.watch(r.next()));
			}
		progress.finish();
		w.close();
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
	 	
	
	public static void main(final String[] args)
		{
		new VcfBurdenFisherH().instanceMainWithExit(args);
		}
	}
