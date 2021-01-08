/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.DoublePredicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

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
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

Variant in that VCF should have one and **only one** ALT allele. Use [https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele](https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele) if needed.

### Output

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
	description="MAF for Cases / Controls ",
	keywords={"vcf","burden","maf","case","control"},
	modificationDate="202000713",
	creationDate="20160418"
	)
public class VcfBurdenMAF
	extends OnePassVcfLauncher
	{
	
	private static final Logger LOG = Logger.build(VcfBurdenMAF.class).make();


	@Parameter(names={"-m","--min-maf","--min-af"},description="select variants where MAF of cases OR MAF of control is greater or equals than min-maf. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double min_AF = 0.0 ;

	@Parameter(names={"-M","--max-maf","--max-af"},description="select variants where MAF of cases OR MAF of control is lower or equal than max-maf. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double max_AF = 0.05 ;
	
	@Parameter(names={"-c","-hr","--hr","--homref"},description="Treat No Call './.' genotypes as HOM_REF '0/0' ")
	private boolean noCallAreHomRef = false;
	
	@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="Don't try to calculate things why variants already FILTERed (faster)")
	private boolean ignoreFiltered=false;
	
	@Parameter(names={"-p","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile=null;
	
	@Parameter(names={"-gtf","--gtf","--gtFiltered"},description="[20180117] Ignore FILTERed **Genotype**")
	private boolean ignore_filtered_genotype=false;
	
	@Parameter(names={"-pfx","--prefix"},description="Prefix for FILTER/INFO. If it is empty and the variant is FILTERed, the variant won'be written to output.")
	private String prefix="Burden";

	
	public VcfBurdenMAF()
		{
		}
	 	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter out)
		{
		final int CASE_POP=0;

		final VCFHeader header0 = in.getHeader();
		
		final String maf_label = "("+this.min_AF+"<= maf <= "+this.max_AF+")";
		final VCFInfoHeaderLine mafCasInfoHeader = (StringUtils.isBlank(this.prefix)?null:new VCFInfoHeaderLine(
				this.prefix + "AF_Cases",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"MAF Cases"
				));
		final VCFInfoHeaderLine mafControlsInfoHeader =  (StringUtils.isBlank(this.prefix)?null:new VCFInfoHeaderLine(
				this.prefix + "AF_Controls",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"MAF Controls"
				));
		
		final VCFInfoHeaderLine acCasInfoHeader =  (StringUtils.isBlank(this.prefix)?null:new VCFInfoHeaderLine(
				this.prefix + "AC_Cases",VCFHeaderLineCount.A,VCFHeaderLineType.Integer,"AC Cases"
				));
		final VCFInfoHeaderLine acControlsInfoHeader =  (StringUtils.isBlank(this.prefix)?null:new VCFInfoHeaderLine(
				this.prefix + "AC_Controls",VCFHeaderLineCount.A,VCFHeaderLineType.Integer,"AC Controls"
				));
		
		final VCFFilterHeaderLine filterCasHeader =  (StringUtils.isBlank(this.prefix)?null:new VCFFilterHeaderLine(
				mafCasInfoHeader.getID(),"MAF of case failed: "+maf_label
				));
		final VCFFilterHeaderLine filterControlsHeader = (StringUtils.isBlank(this.prefix)?null:new VCFFilterHeaderLine(
				mafControlsInfoHeader.getID(),"MAF of controls failed: "+maf_label
				));
		final VCFFilterHeaderLine filterCaseOrControlsHeader =  (StringUtils.isBlank(this.prefix)?null:new VCFFilterHeaderLine(
				this.prefix + "MAFCaseOrControls","MAF of cases OR MAF of controls failed: "+maf_label
				));
	
		final Set<Sample> persons;
		
			{
			try {
				persons = new PedigreeParser().
						parse(this.pedigreeFile).
						getSamplesInVcfHeader(header0).
						filter(S->S.isStatusSet()).
						collect(Collectors.toSet())
						;
				}
			catch(final IOException err)
				{
				LOG.error(err);
				return -1;
				}
			}
			
		final DoublePredicate isInAfRange = AF-> this.min_AF <= AF && AF <= this.max_AF;	
			
		final Set<Sample> caseSamples = persons.stream().
				filter(I->I.isAffected()).
				collect(Collectors.toSet());
		final Set<Sample> controlSamples = persons.stream().
				filter(I->I.isUnaffected()).
				collect(Collectors.toSet());

		if(caseSamples.isEmpty()) {
			LOG.warn("NO case in "+this.pedigreeFile);
		}
		if(controlSamples.isEmpty()) {
			LOG.warn("NO control in "+this.pedigreeFile);
		}
		
		final VCFHeader h2= new VCFHeader(header0);
		
		if(!StringUtils.isBlank(this.prefix)) {
			h2.addMetaDataLine(mafCasInfoHeader);
			h2.addMetaDataLine(mafControlsInfoHeader);
			h2.addMetaDataLine(acCasInfoHeader);
			h2.addMetaDataLine(acControlsInfoHeader);
			h2.addMetaDataLine(filterCasHeader);
			h2.addMetaDataLine(filterControlsHeader);
			h2.addMetaDataLine(filterCaseOrControlsHeader);
			}
		JVarkitVersion.getInstance().addMetaData(this, h2);
		
		out.writeHeader(h2);
		while(in.hasNext())
			{
			final VariantContext ctx = in.next();
			
			if(this.ignoreFiltered && ctx.isFiltered())
				{
				if(!StringUtils.isBlank(this.prefix)) out.add(ctx);
				continue;
				}
			if(!ctx.hasGenotypes())
				{
				if(!StringUtils.isBlank(this.prefix)) out.add(ctx);
				continue;
				}
			
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
					for(final Sample p:(pop==CASE_POP?caseSamples:controlSamples)) 
						{
						/* get genotype for this individual */
						final Genotype genotype = ctx.getGenotype(p.getId());
						if(this.ignore_filtered_genotype && genotype.isFiltered()) continue;
						
						mafCalculator.add(genotype, p.isMale());
						
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
							if(isInAfRange.test(maf))  set_max_maf_cas=false;
							}
						else
							{
							/* add INFO attribute */
							mafCtrlList.add(maf);
							acCtrlList.add((int)ac);
							/* remove FILTER if needed */
							if(isInAfRange.test(maf))  set_max_maf_control=false;
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
			
			if(StringUtils.isBlank(this.prefix)) {
				if(!seen_data || set_max_maf_cas || set_max_maf_control) continue;
				out.add(ctx);
				}
			else
				{
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
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
			}
		return 0;
		}

	@Override
	protected int beforeVcf() {
		if(this.min_AF>this.max_AF) {
			LOG.error("bad values for min/max af");
			return -1;
			}
		return 0;
		}
	

	public static void main(final String[] args)
		{
		new VcfBurdenMAF().instanceMainWithExit(args);
		}
	}
