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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/*
BEGIN_DOC

## Motivation

August 2018, for @MKarakachoff

part of the code was inspired from GATK public code :  https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/annotator/AlleleBalance.java

## Example

```
$ java -jar dist/vcfallelebalance.jar  src/test/resources/test_vcf01.vcf
$ java -jar dist/vcfallelebalance.jar -p src/test/resources/test_vcf01.ped src/test/resources/test_vcf01.vcf
```

END_DOC
 */
@Program(name="vcfallelebalance",
	description="Insert missing allele balance annotation using FORMAT:AD",
	keywords= {"vcf","allele-balance","depth"}
	)
public class VcfAlleleBalance extends Launcher {
	private static final Logger LOG = Logger.build(VcfAlleleBalance.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null; 
	@Parameter(names={"-f","--filtered"},description="ignore FILTER-ed **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"-s","--snp"},description="consider only snps")
	private boolean only_snps = false;
	@Parameter(names={"-p","-ped","--pedigree","--ped"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;

	private Pedigree pedigree;
	private static final String ALLELE_BALANCE_HET_KEY = "ABHet";
	private static final String ALLELE_BALANCE_HOM_KEY = "ABHom";		
	private static final String NON_DIPLOID_RATIO_KEY = "OND";
	private static final String CASE_PREFIX = "CASE_";
	private static final String CTRL_PREFIX = "CTRL_";

	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iterin,
		final VariantContextWriter out
		) {

		final VCFHeader header = iterin.getHeader();
		
		if(this.pedigreeFile!=null)
			{
			try {
				this.pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			}
		else
			{
			this.pedigree = new Pedigree.Parser().parse(header);
			}
		final Set<String> cases_samples = this.pedigree.getAffected().stream().map(P->P.getId()).collect(Collectors.toSet());
		final Set<String> ctr_samples = this.pedigree.getUnaffected().stream().map(P->P.getId()).collect(Collectors.toSet());
		boolean use_dp4_if_ad_missing = false;
		final VCFHeader header2 = new VCFHeader(header);
		if(header.hasGenotypingData()) {
			
			if(header.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)==null) {
				LOG.error("header is issing  FORMAT/"+VCFConstants.GENOTYPE_ALLELE_DEPTHS);
				header2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
				final VCFFormatHeaderLine dp4h = header.getFormatHeaderLine("DP4");
				if(dp4h!=null && dp4h.isFixedCount()  && dp4h.getCount()==4 && dp4h.getType()==VCFHeaderLineType.Integer)
					{
					LOG.warn("I will use FORMAT/DP4 instead of FORMAT/AD");
					use_dp4_if_ad_missing = true;
					}
				}
			
			for(int i=0;i< (this.pedigree.isEmpty()?1:3);i++)
				{
				final String prefix;
				switch(i)
					{
					case 0: prefix="";break;
					case 1: prefix=CASE_PREFIX;break;
					default: prefix=CTRL_PREFIX;break;
					}
				
				String key= prefix+ALLELE_BALANCE_HET_KEY;
				if(header.getInfoHeaderLine(key)!=null)
					{
					LOG.warn("header already contains INFO/"+key);
					}
				else
					{
					header2.addMetaDataLine(
						new VCFInfoHeaderLine(
							key,
							1,
							VCFHeaderLineType.Float,
							prefix.replace("_", ":") + "Allele Balance for heterozygous calls (ref/(ref+alt))"));
					}
				
				key= prefix+ALLELE_BALANCE_HOM_KEY;
				
				if(header.getInfoHeaderLine(key)!=null)
					{
					LOG.warn("header already contains INFO/"+key);
					}
				else
					{
					header2.addMetaDataLine(
						new VCFInfoHeaderLine(
							key,
							1,
							VCFHeaderLineType.Float,
							prefix.replace("_", ":") + "Allele Balance for homozygous calls (A/(A+O)) where A is the allele (ref or alt) and O is anything other"
							));
					}
				
				key= prefix+NON_DIPLOID_RATIO_KEY;
				if(header.getInfoHeaderLine(key)!=null)
					{
					LOG.warn("header already contains INFO/"+key);
					}
				else
					{
					header2.addMetaDataLine(
						new VCFInfoHeaderLine(
							key,
							1,
							VCFHeaderLineType.Float,
							prefix.replace("_", ":") + "Overall non-diploid ratio (alleles/(alleles+non-alleles))"
							));
					}
				}
			}
		JVarkitVersion.getInstance().addMetaData(VcfAlleleBalance.class.getSimpleName(), header2);
		out.writeHeader(header2);
		
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.
				newInstance().
				logger(LOG).
				dictionary(header).
				build();
		while(iterin.hasNext())
			{
			VariantContext ctx = progress.apply(iterin.next());
			if(!(ctx.hasGenotypes() || ctx.isBiallelic()))
				{
				out.add(ctx);
				continue;
				}
			
			if(use_dp4_if_ad_missing && ctx.getNAlleles()==2)
				{
				final List<Genotype> newgt = new ArrayList<>(ctx.getNSamples());
				for(final Genotype gt: ctx.getGenotypes())
					{
					if(!gt.hasAnyAttribute("DP4") || gt.hasAD() || gt.isHetNonRef())
						{
						newgt.add(gt);
						}
					else
						{
						Object o= gt.getAnyAttribute("DP4");
						if(o==null || !(o instanceof List) || List.class.cast(o).size()!=4) continue;
						final List<?> dp4 = (List<?>)o;
						final int ad[] = new int[ctx.getNAlleles()];
						Arrays.fill(ad, 0);
						
						//Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases
						ad[0] = Integer.class.cast(dp4.get(0)) + Integer.class.cast(dp4.get(1));
						ad[1] = Integer.class.cast(dp4.get(2)) + Integer.class.cast(dp4.get(3));
						
						newgt.add(new GenotypeBuilder(gt).AD(ad).make());
						}
					}
				
				ctx = new VariantContextBuilder(ctx).genotypes(newgt).make();
				}
			
			
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			computeAB(ctx,vcb,(GT)->true,"");
			if(!pedigree.isEmpty())
				{
				computeAB(ctx,vcb,
						(GT)->cases_samples.contains(GT.getSampleName()),
						CASE_PREFIX
						);
				computeAB(ctx,vcb,
						(GT)->ctr_samples.contains(GT.getSampleName()),
						CTRL_PREFIX
						);
				}
			out.add(vcb.make());
			}
		progress.close();
		return 0;
		}
	
	private VariantContextBuilder computeAB(
		final VariantContext ctx,
		final VariantContextBuilder vcb,
		final Predicate<Genotype> acceptGt,
		final String prefix
		) {

        vcb.rmAttribute(prefix + ALLELE_BALANCE_HET_KEY);
        vcb.rmAttribute(prefix + ALLELE_BALANCE_HOM_KEY);
        vcb.rmAttribute(prefix + NON_DIPLOID_RATIO_KEY);

		if(!ctx.hasGenotypes()) return vcb;
        double refCountInHetSamples = 0.0;
        double altCountInHetSamples = 0.0;
        double correctCountInHomSamples = 0.0;
        double incorrectCountInHomSamples = 0.0;
        double nonDiploidCount = 0.0;
        double totalReadCount = 0.0;

        for ( final Genotype genotype : ctx.getGenotypes() ) {
        	if(!acceptGt.test(genotype)) continue;
        	if(this.only_snps && !ctx.isSNP()) continue;
        	if(genotype.isNoCall()) continue;
        	if(!genotype.hasAD()) continue;
        	if(this.ignore_filtered_genotypes && genotype.isFiltered()) continue;
            final int[] alleleCounts = genotype.getAD();
            if (alleleCounts == null) continue;
            if(alleleCounts.length!=ctx.getNAlleles()) {
            	LOG.warn("expected AD.length="+ctx.getNAlleles()+" for "+genotype.getSampleName()+" at "+ctx.getContig()+":"+ctx.getStart());
            	continue;
            }
            
            // inspired from GATK/Broad below :  https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/annotator/AlleleBalance.java
            
            final long totalReads = IntStream.of(alleleCounts).mapToLong(I->I).sum();
            
            if ( genotype.isHet() ) {
                refCountInHetSamples += alleleCounts[0];
                altCountInHetSamples += alleleCounts[1];
                nonDiploidCount += totalReads - (alleleCounts[0] + alleleCounts[1]);
                totalReadCount += totalReads;
            	} 
            else if ( genotype.isHom() ) {
                final int alleleIndex = genotype.isHomRef() ?  0 : 1 ;
                final int alleleCount = alleleCounts[alleleIndex];
                int bestOtherCount = 0;
                for(int n = 0; n < alleleCounts.length; n++){
                    if( n != alleleIndex && alleleCounts[n] > bestOtherCount ) {
                        bestOtherCount = alleleCounts[n];
                    }
                }
                correctCountInHomSamples += alleleCount;
                incorrectCountInHomSamples += bestOtherCount;
                nonDiploidCount += totalReads - alleleCount;
                totalReadCount += totalReads;
            }
        }
    final double diploidCountInHetSamples = altCountInHetSamples + refCountInHetSamples;
    final double diploidCountInHomSamples = correctCountInHomSamples + incorrectCountInHomSamples;
    
    if ( diploidCountInHetSamples > 0.0 ) {
    	 vcb.attribute(prefix + ALLELE_BALANCE_HET_KEY, refCountInHetSamples / diploidCountInHetSamples);
    }

    if ( diploidCountInHomSamples > 0.0 ) {
    	vcb.attribute(prefix + ALLELE_BALANCE_HOM_KEY,  correctCountInHomSamples / diploidCountInHomSamples);
    }

    if ( totalReadCount > 0.0 ) {
    	vcb.attribute(prefix + NON_DIPLOID_RATIO_KEY, nonDiploidCount / totalReadCount);
    	}
    return vcb;

	}
	
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args, this.outputFile);
		}
	public static void main(final String[] args) {
		new VcfAlleleBalance().instanceMainWithExit(args);
	}

}
