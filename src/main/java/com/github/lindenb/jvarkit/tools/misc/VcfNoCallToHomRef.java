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


History:
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

/**
 BEGIN_DOC

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
		biostars={276811}
		)
public class VcfNoCallToHomRef extends  Launcher
	{
	private static final Logger LOG=Logger.build(VcfNoCallToHomRef.class).make();
	
	@Parameter(names={"-o","--out"},required=false,description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output=null;

	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();
	
	public VcfNoCallToHomRef()
		{
		}
	
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
		{
		@Parameter(names={"-s","--includeSamples"},description="only converts those samples. Default: all samples are converted.")
		private  Set<String> includeSamples = new HashSet<>();
		@Parameter(names={"-sf","--includeSamplesFile"},description="only converts those samples. Default: all samples are converted. One sample per line.")
		private File includeSampleFile = null;
		@Parameter(names={"-dp","--depth"},description="Default DEPTH. negative = don't set depth.")
		private int defaultDepth=10;
		@Parameter(names={"-gq","--gq","--GT"},description="Default Genotype quality: negative : don't set GQ.")
		private int defaultGT=1;
		@Parameter(names={"-f","--filter"},description="Set this Genotype FILTER for converted genotype")
		private String filterName = null;
		@Parameter(names={"-p","--ploidy"},description="ploidy")
		private int ploidy = 2;
		@Parameter(names={"--noRecount"},description="do not recount DP/AC/AN/AF atttributes")
		private boolean doNotRecountAcAnAfAttributes = false;
		@Parameter(names={"-x"},description="do not recount DP/AC/AN/AF atttributes")
		private boolean doNotFix = false;

		
		private final Set<String> allIncludeSamples = new HashSet<>();
		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			/** excluded samples */
			private final Set<String> included = new HashSet<>(CtxWriterFactory.this.allIncludeSamples);
			private long countFixedGenotypes = 0L;
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				}
			@Override
			public void writeHeader(final VCFHeader header) {
				final VCFHeader header2 = new VCFHeader(header);
				if(this.included.isEmpty())
					{
					this.included.addAll(header.getSampleNamesInOrder());
					}
				else
					{
					this.included.retainAll(header.getSampleNamesInOrder());
					}
				if(!StringUtil.isBlank(CtxWriterFactory.this.filterName))
					{
					final VCFFilterHeaderLine fh = new VCFFilterHeaderLine(
							CtxWriterFactory.this.filterName,
							"NOCALL Genotypes converted to HOM_REF by "+VcfNoCallToHomRef.class.getName()
							);
					header2.addMetaDataLine(fh);
					}
		
				for(final String key: new String[]{VCFConstants.DEPTH_KEY,VCFConstants.ALLELE_COUNT_KEY,VCFConstants.ALLELE_FREQUENCY_KEY,VCFConstants.ALLELE_NUMBER_KEY})
					{
					if(!header.hasInfoLine(key))
						{
						header2.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(key));
						}
					}
				for(final String key: new String[]{VCFConstants.GENOTYPE_FILTER_KEY,VCFConstants.DEPTH_KEY,VCFConstants.GENOTYPE_KEY,VCFConstants.GENOTYPE_QUALITY_KEY})
					{
					if(!header.hasFormatLine(key))
						{
						header2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(key));
						}
					}
				super.writeHeader(header2);
				}
			@Override
			public void add(final VariantContext ctx) {
			
				final Map<String,Genotype> sample2genotypes = new HashMap<>(ctx.getNSamples());
				
				/** hom-ref alleles for GT builder */
				final List<Allele> homRefAlleles = new ArrayList<>(CtxWriterFactory.this.ploidy);
				for(int i=0;i< CtxWriterFactory.this.ploidy;++i)
					{
					homRefAlleles.add(ctx.getReference());
					}

				
				ctx.getGenotypes().stream().forEach(G->{
					if(G.isCalled() || !this.included.contains(G.getSampleName()))
						{
						sample2genotypes.put(G.getSampleName(), G);
						return;
						}

					final GenotypeBuilder gb = new GenotypeBuilder(G.getSampleName(), homRefAlleles);
					if( CtxWriterFactory.this.defaultDepth >= 0)
						{
						gb.DP(CtxWriterFactory.this.defaultDepth );
						}
					if( CtxWriterFactory.this.defaultGT >= 0)
						{
						gb.GQ(CtxWriterFactory.this.defaultGT );
						}
					if(!StringUtil.isBlank(CtxWriterFactory.this.filterName))
						{
						gb.filter(CtxWriterFactory.this.filterName);
						}
					countFixedGenotypes++;
					sample2genotypes.put(G.getSampleName(), gb.make());
					});
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.genotypes(sample2genotypes.values());
				if(!doNotRecountAcAnAfAttributes)
					{
					vcb.attribute(VCFConstants.DEPTH_KEY, sample2genotypes.values().stream().filter(G->G.hasDP()).mapToInt(G->G.getDP()).sum());
					final int an = (int)sample2genotypes.values().stream().flatMap(G->G.getAlleles().stream()).filter(A->A.isCalled()).count();
					vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,an);
					final List<Integer> ac = new ArrayList<>();
					final List<Float> af = new ArrayList<>();
					for(final Allele alt: ctx.getAlternateAlleles())
						{
						final int acn = (int)sample2genotypes.values().stream().
								flatMap(G->G.getAlleles().stream()).
								filter(A->A.equals(alt)).
								count();
						ac.add(acn);
						af.add(an<=0?null:new Float((float)acn/(float)an));
						}
					if(!ac.isEmpty()) {
						vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,ac);
						vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,af);
					} else
						{
						vcb.rmAttribute(VCFConstants.ALLELE_COUNT_KEY);
						vcb.rmAttribute(VCFConstants.ALLELE_FREQUENCY_KEY);
						}
					}
				
				super.add(vcb.make());
				}
			@Override
			public void close() {
				LOG.info("Number of fixed genotypes : "+countFixedGenotypes);
				super.close();
				}
			}
		
		@Override
		public int initialize() {
			if(this.ploidy<1) {
				LOG.error("bad ploidy <1");
				return -1;
			}
			this.allIncludeSamples.addAll(this.includeSamples);
			if(this.includeSampleFile!=null) {
				try {
					this.allIncludeSamples.addAll(
						Files.lines(this.includeSampleFile.toPath()).
							filter(L->!(L.isEmpty()|| L.startsWith("#"))).
							collect(Collectors.toSet())
						);
					}
				catch(final IOException err)
					{
					LOG.error(err);
					return -1;
					}
				}
			return 0;
			}
		
		@Override
		public CtxWriter open(final VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}
		@Override
		public void close() throws IOException {
			this.allIncludeSamples.clear();
			}
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter delegate
			) 
		{
		final CtxWriterFactory.CtxWriter out = this.component.open(delegate);
		out.writeHeader(in.getHeader());
		final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
		while(in.hasNext())
			{
			out.add(progress.watch(in.next()));
			}
		progress.finish();
		out.close();
		return 0;
		}
	
	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		return new PostponedVariantContextWriter(this.writingVcfArgs,stdout(),outorNull);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.output);
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
		
	public static void main(final String[] args)
		{
		new VcfNoCallToHomRef().instanceMainWithExit(args);
		}
	}
