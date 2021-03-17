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


*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;


import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import htsjdk.variant.vcf.VCFIterator;

/**
 BEGIN_DOC

## Deprecated


deprecated. Use `bcftools +setGT`

## Example

original VCF

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	.	GT:PL	./.	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	.	GT:PL	0/0:0,255,133	0/1:40,0,31	./.	./.
```

default invocation

```
$ java -jar dist/vcfnocall2homref.jar input.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC=2;AF=0.25;AN=8;DP=10	GT:DP:GQ:PL	0/0:10:1	0/0:.:.:0,255,127	0/0:.:.:0,255,137	1/1:.:.:70,255,0
rotavirus	91	.	A	T	5.45	.	AC=1;AF=0.125;AN=8;DP=20	GT:DP:GQ:PL	0/0:.:.:0,255,133	0/1:.:.:40,0,310/0:10:1	0/0:10:1
```

convert S3 and S4 only

```
$ java -jar dist/vcfnocall2homref.jar  -f CONVERTED -s S3 -s S4  ~/jeter.vcf 
##fileformat=VCFv4.2
##FILTER=<ID=CONVERTED,Description="NOCALL Genotypes converted to HOM_REF by com.github.lindenb.jvarkit.tools.misc.VcfNoCallToHomRef">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC=2;AF=0.33333334;AN=6;DP=0	GT:PL	./.	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC=1;AF=0.125;AN=8;DP=20	GT:DP:FT:GQ:PL	0/0:.:PASS:.:0,255,133	0/1:.:PASS:.:40,0,31	0/0:10:CONVERTED:1	0/0:10:CONVERTED:1
```

END_DOC
 */
@Program(
	name="vcfnocall2homref",
	description="Convert the UNCALLED gentoypes in a VCF to HOM_REF. This tool can be used after using GATK CombineVariants.",
	keywords={"vcf"},
	creationDate="20170914",
	modificationDate="20200720",
	biostars={276811},
	deprecatedMsg="use bcftools plugin: +setGT"
	)
public class VcfNoCallToHomRef extends OnePassVcfLauncher
	{
	private static final Logger LOG=Logger.build(VcfNoCallToHomRef.class).make();
		
	@Parameter(names={"-s","--includeSamples"},description="only converts those samples. Default: all samples are converted.")
	private  Set<String> includeSamples = new HashSet<>();
	@Parameter(names={"-sf","--includeSamplesFile"},description="only converts those samples. Default: all samples are converted. One sample per line.")
	private Path includeSampleFile = null;
	@Parameter(names={"-dp","--depth"},description="Default DEPTH. negative = don't set depth.")
	private int defaultDepth=10;
	@Parameter(names={"-gq","--gq","--GQ"},description="Default Genotype quality: negative : don't set GQ.")
	private int defaultGQ=1;
	@Parameter(names={"-f","--filter"},description="Set this **Genotype** FILTER for converted genotype")
	private String filterName = null;
	@Parameter(names={"-p","--ploidy"},description="ploidy")
	private int ploidy = 2;	
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter out
			) 
		{
		final Set<String> allIncludeSamples = new HashSet<>(this.includeSamples);
		
		if(this.includeSampleFile!=null) {
			try {
				allIncludeSamples.addAll(
					Files.lines(this.includeSampleFile).
						filter(L->!(StringUtils.isBlank(L)|| L.startsWith("#"))).
						collect(Collectors.toSet())
					);
				}
			catch(final IOException err)
				{
				LOG.error(err);
				return -1;
				}
			}
		
		final VCFHeader header = in.getHeader();
		final VCFHeader header2 = new VCFHeader(header);
		final Set<VCFHeaderLine> headerLines = new HashSet<>();
		for(final String key : new String[] {
				VCFConstants.GENOTYPE_ALLELE_DEPTHS,
				VCFConstants.GENOTYPE_FILTER_KEY,
				VCFConstants.GENOTYPE_QUALITY_KEY,
				VCFConstants.DEPTH_KEY
				})
		if(!header2.hasFilterLine(key)) {
			VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,key);
			}
		
		headerLines.forEach(H->header2.addMetaDataLine(H));

		final VariantAttributesRecalculator recalc = new VariantAttributesRecalculator();
		recalc.setHeader(header2);
		
		
		
		// change all
		if(this.includeSamples.isEmpty() && this.includeSampleFile==null)
			{
			allIncludeSamples.addAll(header.getSampleNamesInOrder());
			}
		// only those in the vcf
		else
			{
			allIncludeSamples.retainAll(header.getSampleNamesInOrder());
			}

		JVarkitVersion.getInstance().addMetaData(this, header2);
		out.writeHeader(header2);
		long countFixedGenotypes = 0L;
		while(in.hasNext())
			{
			final VariantContext ctx = in.next();				

			final List<Genotype> sample2genotypes = new ArrayList<>(ctx.getNSamples());
			
			/** hom-ref alleles for GT builder */
			final List<Allele> homRefAlleles = new ArrayList<>(this.ploidy);
			for(int i=0;i< this.ploidy;++i)
				{
				homRefAlleles.add(ctx.getReference());
				}

			
			for(final Genotype G:ctx.getGenotypes()) {
				if(G.isCalled() || !allIncludeSamples.contains(G.getSampleName()))
					{
					sample2genotypes.add(G);
					}
				else
					{
					final GenotypeBuilder gb = new GenotypeBuilder(G);
					gb.alleles(homRefAlleles);
					if( this.defaultDepth >= 0)
						{
						gb.DP(this.defaultDepth );
						int ad[]=new int[ctx.getNAlleles()];
						Arrays.fill(ad, 0);
						ad[0] = this.defaultDepth;
						gb.AD(ad);
						}
					if( this.defaultGQ >= 0)
						{
						gb.GQ(this.defaultGQ );
						}
					if(!StringUtil.isBlank(this.filterName))
						{
						gb.filter(this.filterName);
						}
					countFixedGenotypes++;
					sample2genotypes.add(gb.make());
					}
				};
			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.genotypes(sample2genotypes);
			out.add(recalc.apply(vcb.make()));
			}
		LOG.info("Number of fixed genotypes : "+countFixedGenotypes);
		return 0;
		}

		
	public static void main(final String[] args)
		{
		new VcfNoCallToHomRef().instanceMainWithExit(args);
		}
	}
