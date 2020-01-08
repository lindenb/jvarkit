/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.lumpysv.LumpyConstants;
import com.github.lindenb.jvarkit.util.Pedigree;
import htsjdk.variant.vcf.VCFIterator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

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

# Example

```
bcftools annotate -x 'INFO' src/test/resources/rotavirus_rf.ann.vcf.gz |\
	java -jar dist/vcfburdenmaf.jar --pedigree pedigree.ped

##fileformat=VCFv4.2
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##FILTER=<ID=BurdenMAFCas,Description="MAF of cases is greater than 0.05">
##FILTER=<ID=BurdenMAFCaseOrControls,Description="MAF of (cases OR controls) is greater than 0.05">
##FILTER=<ID=BurdenMAFControls,Description="MAF of controls is greater than 0.05">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=BurdenMAFCas,Number=A,Type=Float,Description="Burden Filter F2. MAF Cases">
##INFO=<ID=BurdenMAFControls,Number=A,Type=Float,Description="Burden Filter F2. MAF Controls">
##contig=<ID=RF01,length=3302>
##contig=<ID=RF02,length=2687>
##contig=<ID=RF03,length=2592>
##contig=<ID=RF04,length=2362>
##contig=<ID=RF05,length=1579>
##contig=<ID=RF06,length=1356>
##contig=<ID=RF07,length=1074>
##contig=<ID=RF08,length=1059>
##contig=<ID=RF09,length=1062>
##contig=<ID=RF10,length=751>
##contig=<ID=RF11,length=666>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	970	.	A	C	48.67	.	BurdenMAFCas=0.00;BurdenMAFControls=0.00	GT:PL	0/0:0,9,47	0/0:0,18,73	0/0:0,18,73	0/0:0,33,116	1/1:95,24,0
RF02	251	.	A	T	21.29	BurdenMAFCaseOrControls;BurdenMAFControls	BurdenMAFCas=0.00;BurdenMAFControls=0.500	GT:PL	0/0:0,15,57	0/1:31,0,5	0/1:31,0,5	0/0:0,9,42	0/0:0,24,69
RF02	578	.	G	A	53	.	BurdenMAFCas=0.00;BurdenMAFControls=0.00	GT:PL	0/0:0,33,122	0/0:0,39,135	0/0:0,39,135	1/1:100,30,0	0/0:0,27,109
RF02	877	.	T	A	3.45	BurdenMAFCas;BurdenMAFCaseOrControls	BurdenMAFCas=0.500;BurdenMAFControls=0.00	GT:PL	0/1:37,0,50	0/0:0,22,116	0/0:0,22,116	0/0:0,21,94	0/0:0,12,62
RF02	1726	.	T	G	8.23	BurdenMAFCaseOrControls;BurdenMAFControls	BurdenMAFCas=0.00;BurdenMAFControls=0.500	GT:PL	0/0:0,18,83	0/1:24,0,40	0/1:24,0,40	0/0:0,27,111	0/0:0,10,78
RF02	1962	.	TACA	TA	33.43	BurdenMAFCas;BurdenMAFCaseOrControls	BurdenMAFCas=0.500;BurdenMAFControls=0.00	GT:PL	0/1:70,0,159	0/0:0,15,225	0/0:0,15,225	0/0:0,27,231	0/0:0,27,168
```

END_DOC

*/
@Program(name="vcfburdenmaf",
	description="Burden : MAF for Cases / Controls ",
	keywords={"vcf","burden","maf","case","control"},
	modificationDate="20190628"
	)
public class VcfBurdenMAF
	extends Launcher
	{
	private static final int CASE_POP=0;
	
	private static final Logger LOG = Logger.build(VcfBurdenMAF.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-maxMAF","--maxMAF","-maxAF","--maxAF"},description="if MAF of cases OR MAF of control is greater than maxMAF, the the FILTER Column is Filled")
	private double maxMAF = 0.05 ;
	
	@Parameter(names={"-c","--homref"},description="Treat No Call './.' genotypes as HomRef")
	private boolean noCallAreHomRef = false;
	
	@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Don't try to calculate things why variants already FILTERed (faster)")
	private boolean ignoreFiltered=false;
	
	@Parameter(names={"-p","--pedigree"},description="[20180117] Pedigree file. Default: use the pedigree data in the VCF header." + Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile=null;
	
	@Parameter(names={"-gtf","--gtf","--gtFiltered"},description="[20180117] Ignore FILTERed **Genotype**")
	private boolean ignore_filtered_genotype=false;
	
	@Parameter(names={"-lumpy-su-min","--lumpy-su-min"},description="[20180117] if variant identified as LUMPy-SV variant. This is the minimal number of 'SU' to consider the genotype as a variant.")
	private int lumpy_SU_threshold=1;

	@Parameter(names={"-pfx","--prefix"},description="Prefix for FILTER/INFO")
	private String prefix="Burden";

	
	public VcfBurdenMAF()
		{
		}
	 	
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter out)
		{
		final VCFHeader header0 = in.getHeader();
		final ProgressFactory.Watcher<VariantContext> progess = ProgressFactory.newInstance().dictionary(header0).logger(LOG).build();
		
		final Function<VCFHeader,Set<Pedigree.Person>> caseControlExtractor = 
				(header)->  new Pedigree.CaseControlExtractor().extract(header0);


		final VCFInfoHeaderLine mafCasInfoHeader = new VCFInfoHeaderLine(
				this.prefix + "AF_Cases",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"MAF Cases"
				);
		final VCFInfoHeaderLine mafControlsInfoHeader = new VCFInfoHeaderLine(
				this.prefix + "AF_Controls",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"MAF Controls"
				);
		
		final VCFInfoHeaderLine acCasInfoHeader = new VCFInfoHeaderLine(
				this.prefix + "AC_Cases",VCFHeaderLineCount.A,VCFHeaderLineType.Integer,"AC Cases"
				);
		final VCFInfoHeaderLine acControlsInfoHeader = new VCFInfoHeaderLine(
				this.prefix + "AC_Controls",VCFHeaderLineCount.A,VCFHeaderLineType.Integer,"AC Controls"
				);
		
		final VCFFilterHeaderLine filterCasHeader = new VCFFilterHeaderLine(
				mafCasInfoHeader.getID(),"MAF of cases is greater than "+this.maxMAF
				);
		final VCFFilterHeaderLine filterControlsHeader = new VCFFilterHeaderLine(
				mafControlsInfoHeader.getID(),"MAF of controls is greater than "+this.maxMAF
				);
		final VCFFilterHeaderLine filterCaseOrControlsHeader = new VCFFilterHeaderLine(
				this.prefix + "MAFCaseOrControls","MAF of (cases OR controls) is greater than "+this.maxMAF
				);			
		

		
		final boolean is_lumpy_vcf_header = LumpyConstants.isLumpyHeader(header0);
		final Set<Pedigree.Person> persons;
		if( this.pedigreeFile == null)
			{
			persons= caseControlExtractor.apply(header0);
			}
		else
			{
			try {
				persons = new Pedigree.CaseControlExtractor().extract(
						header0,
						new Pedigree.Parser().parse(this.pedigreeFile)
						);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			}
		final Set<Pedigree.Person> caseSamples = persons.stream().
				filter(I->I.isAffected()).collect(Collectors.toSet());
		final Set<Pedigree.Person> controlSamples = persons.stream().
				filter(I->I.isUnaffected()).collect(Collectors.toSet());

		final VCFHeader h2= new VCFHeader(header0);
		h2.addMetaDataLine(mafCasInfoHeader);
		h2.addMetaDataLine(mafControlsInfoHeader);
		h2.addMetaDataLine(acCasInfoHeader);
		h2.addMetaDataLine(acControlsInfoHeader);
		h2.addMetaDataLine(filterCasHeader);
		h2.addMetaDataLine(filterControlsHeader);
		h2.addMetaDataLine(filterCaseOrControlsHeader);
		
		out.writeHeader(h2);
		while(in.hasNext())
			{
			final VariantContext ctx = progess.apply(in.next());
			
			if(this.ignoreFiltered && ctx.isFiltered())
				{
				out.add(ctx);
				continue;
				}
			if(!ctx.hasGenotypes())
				{
				out.add(ctx);
				continue;
				}
			
			final boolean identified_as_lumpy= 
					is_lumpy_vcf_header && 
					LumpyConstants.isLumpyVariant(ctx)
					;
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			final List<Double> mafCasList = new ArrayList<>(); 
			final List<Double> mafCtrlList = new ArrayList<>(); 
			final List<Integer> acCasList = new ArrayList<>(); 
			final List<Integer> acCtrlList = new ArrayList<>(); 
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
						if(this.ignore_filtered_genotype && genotype.isFiltered()) continue;
						
						/* this is a lumpy genotype */
						if(identified_as_lumpy)
							{
							if(!genotype.hasExtendedAttribute("SU"))
								{
								throw new JvarkitException.FileFormatError(
										"Variant identified as lumpysv, but not attribute 'SU' defined in genotye "+genotype);
								}
							@SuppressWarnings("deprecation")
							final int su_count = genotype.getAttributeAsInt("SU", 0);
							final boolean genotype_contains_allele = su_count>= this.lumpy_SU_threshold;
							mafCalculator.add(
									new GenotypeBuilder(genotype.getSampleName(),
											genotype_contains_allele ?
											Arrays.asList(observed_alt,observed_alt):
											Arrays.asList(ctx.getReference(),ctx.getReference())
											).make()
									, p.isMale());
							}
						else /* this is a not a lumpy genotype , regular case...*/
							{
							mafCalculator.add(genotype, p.isMale());
							}
						/* if(pop==CASE_POP && genotype.isCalled()) LOG.info("DEBUGMAF: "+p+" "+genotype); */
						}/* end of loop over persons */
					/* at least one genotype found */
					if(!mafCalculator.isEmpty())
						{
						seen_data=true;
						
						/* get MAF */
						final double maf = mafCalculator.getMaf();
						final double ac = mafCalculator.getCountAlt();
						if(pop == CASE_POP) {
							/* add INFO attribute */
							mafCasList.add(maf);
							acCasList.add((int)ac);
							/* remove FILTER if needed */
							if(maf<=this.maxMAF)  set_max_maf_cas=false;
							}
						else
							{
							/* add INFO attribute */
							mafCtrlList.add(maf);
							acCtrlList.add((int)ac);
							/* remove FILTER if needed */
							if(maf<=this.maxMAF)  set_max_maf_control=false;
							}
						} 
					else
						{
						if(pop == CASE_POP) {
							mafCasList.add(0.0);
							acCasList.add(0);
							set_max_maf_cas=false;
						} else
						{
							mafCtrlList.add(0.0);
							acCtrlList.add(0);
							set_max_maf_control=false;
						}
						}
					}/* end of loop over pop */
				}/* end loop over alt allele */
			
			
			vcb.attribute(mafCasInfoHeader.getID(),mafCasList);
			vcb.attribute(mafControlsInfoHeader.getID(),mafCtrlList);
			vcb.attribute(acCasInfoHeader.getID(),acCasList);
			vcb.attribute(acControlsInfoHeader.getID(),acCtrlList);
			
			if(seen_data) {
				if(set_max_maf_cas) vcb.filter(filterCasHeader.getID());
				if(set_max_maf_control) vcb.filter(filterControlsHeader.getID());
				if(set_max_maf_cas || set_max_maf_control) {
					vcb.filter(filterCaseOrControlsHeader.getID());
					}
				}

			out.add(vcb.make());
			}
		progess.close();
		return 0;
		}

	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	 	
	
	public static void main(final String[] args)
		{
		new VcfBurdenMAF().instanceMainWithExit(args);
		}
	}
