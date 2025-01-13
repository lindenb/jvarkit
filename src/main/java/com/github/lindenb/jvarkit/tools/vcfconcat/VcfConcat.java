/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Motivation

concat numerous VCF, without checking files at start, doesn't sort.

## Input

input is a set of VCF files or a file with the suffix .list containing the path to the vcfs

## Example

### From stdin

```bash
$ find ./ -name "*.vcf" | grep Sample1 | java -jar dist/jvarkit vcfconcat > out.vcf
```

### From files

```bash
$ java -jar dist/jvarkit vcfconcat Sample1.samtools.vcf Sample1.gatk.vcf > out.vcf
$ java -jar dist/jvarkit vcfconcat paths.list > out.vcf
```


## See also

bcftools concat

END_DOC
*/
@Program(name="vcfconcat",
	keywords={"vcf"},
	creationDate = "20131230",
	modificationDate = "20240426",
	description="Concatenate VCFs with same sample. See also bcftools concat",
	generate_doc = true,
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfConcat extends Launcher
	{
	private static final Logger LOG =Logger.build(VcfConcat.class).make();
	private enum SamplePeek {none,all,with_alt};
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-T","--tag"},description="if not empty, add INFO/tag containing the source/path of the variant")
	private String variantsourceTag="";
	@Parameter(names={"-G","--drop-genotypes"},description="Drop genotypes")
	private boolean drop_genotypes = false;
	@Parameter(names={"-S","--samples"},description="implies --drop-genotypes")
	private SamplePeek samplePeek = SamplePeek.none;
	
	@Parameter(names={"--chrom","--contig"},description="limit to that chromosome")
	private String limitChrom = null;

	@Parameter(names={"--merge"},description="merge all samples. First Scan all files to get all distinct samples")
	private boolean merge_distinct_samples = false;

	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();

	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		try
			{
			final Predicate<Genotype> genotype_peeker;
			if(!this.samplePeek.equals(SamplePeek.none)) {
				this.drop_genotypes = true;
				}
			switch(this.samplePeek) {
				case all : genotype_peeker = G->true; break;
				case with_alt : genotype_peeker = G->G.getAlleles().stream().anyMatch(A->!(A.isReference() || A.isNoCall())); break;
				default : genotype_peeker = G->false; break;
				}
			
			final List<Path> vcfs;
			
			
			if(args.isEmpty())
				{
				try(BufferedReader br=IOUtils.openStreamForBufferedReader(stdin())) {
					vcfs = br.lines().
						filter(line->!(StringUtils.isBlank(line) || line.startsWith("#"))).
						map(F->Paths.get(F)).
						collect(Collectors.toList());
					}
				}
			else
				{
				vcfs = IOUtils.unrollPaths(args);
				}
			
			if(vcfs.isEmpty())
				{
				LOG.error("No input");
				return -1;
				}
			
			final Set<String> distinct_samples;
			if(!drop_genotypes && this.merge_distinct_samples) {
				final Set<String> set = new TreeSet<>();
				for(Path path: vcfs) {
					set.addAll(VcfHeaderExtractor.decode(path).getGenotypeSamples());
					}
				distinct_samples = set.isEmpty()?null:set;
				}
			else
				{
				distinct_samples = null;
				}
			
			VCFInfoHeaderLine variantSourceHeader = null;
			VCFInfoHeaderLine variantSampleHeader = null;
			VCFHeader firstHeader=null;
			long count_variants=0L;
			SAMSequenceDictionary firstDict=null;
			final long initMilliSec = System.currentTimeMillis();
			for(int i=0;i< vcfs.size();i++) {
				final Path vcfPath = vcfs.get(i);
				final long startMilliSec = System.currentTimeMillis();
				LOG.info(String.valueOf(i+1)+"/"+vcfs.size()+" "+vcfPath);
				try(VCFIterator in = new BcfIteratorBuilder().open(vcfPath)) {
					final VCFHeader header0 = in.getHeader();
					final VCFHeader header = drop_genotypes?
							new VCFHeader(header0.getMetaDataInInputOrder()):
							(distinct_samples==null?header0:new VCFHeader(header0.getMetaDataInInputOrder(),distinct_samples));
					final SAMSequenceDictionary dict= header.getSequenceDictionary();
					if(firstHeader==null) {
						w=this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
						firstHeader = header;
						firstDict = dict;
						if(dict!=null && limitChrom!=null && dict.getSequence(this.limitChrom)==null) {
							throw new JvarkitException.ContigNotFoundInDictionary(limitChrom, dict);
							}
						if(!StringUtils.isBlank(this.variantsourceTag)) {
							variantSourceHeader = new VCFInfoHeaderLine(
									this.variantsourceTag,
									1,VCFHeaderLineType.String,
									"Origin File of Variant"
									);
							header.addMetaDataLine(variantSourceHeader);
							}
						
						if(!this.samplePeek.equals(SamplePeek.none)) {
							variantSampleHeader = new VCFInfoHeaderLine(
									"SAMPLES_"+this.samplePeek.name().toUpperCase(),
									VCFHeaderLineCount.UNBOUNDED,
									VCFHeaderLineType.String,
									"Samples (filtered with "+ this.samplePeek.name()+")"
									);
							header.addMetaDataLine(variantSampleHeader);
							}
						
						JVarkitVersion.getInstance().addMetaData(this, header);
						w.writeHeader(firstHeader);
						}
					else
						{
						if(firstDict!=null && dict!=null) SequenceUtil.assertSequenceDictionariesEqual(firstDict, dict);
						if(!this.drop_genotypes && !firstHeader.getGenotypeSamples().equals(header.getGenotypeSamples())) {
							LOG.error("Samples names/order mismatch between "+ vcfs.get(0)+" and "+ vcfPath+
									". You can also use --merge to merge all genotypes, or --drop-genotypes");
							return -1;
							}
						}
					while(in.hasNext()) {
						VariantContext ctx = in.next();
						if(limitChrom!=null && !ctx.getContig().equals(limitChrom)) {
							continue;
							}
												
						if(drop_genotypes || variantSourceHeader!=null) {
							final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
							if(variantSourceHeader!=null) vcb.attribute(
									variantSourceHeader.getID(),
									VCFUtils.escapeInfoField(vcfPath.toString())
									);
							
							if(variantSampleHeader!=null && ctx.hasGenotypes()) {
								final Set<String> samples = ctx.getGenotypes().stream().
											filter(genotype_peeker).
											map(G->G.getSampleName()).
											collect(Collectors.toCollection(TreeSet::new));
								if(!samples.isEmpty()) {
									vcb.attribute(
											variantSampleHeader.getID(),
										new ArrayList<String>(samples)
										);
									}
								}
							if(drop_genotypes) {
								vcb.noGenotypes();
								}
							
							ctx = vcb.make();
							}
						
						if(distinct_samples!=null)
							{
							final int ploidy = ctx.getGenotypes().stream().mapToInt(G->G.getPloidy()).max().orElse(2);
							final List<Genotype> genotypes=new ArrayList<>(distinct_samples.size());
							for(String sn:distinct_samples) {
								Genotype g = ctx.getGenotype(sn);
								if(g==null) {
									g = GenotypeBuilder.createMissing(sn, ploidy);
									}
								genotypes.add(g);
								}
							ctx= new VariantContextBuilder(ctx).genotypes(genotypes).make();
							}
						
						w.add(ctx);
						count_variants++;
						}
					final long millisecPerVcf  = (System.currentTimeMillis() - initMilliSec)/(i+1L);
					LOG.info("N="+count_variants+". That took: "+StringUtils.niceDuration(System.currentTimeMillis() - startMilliSec)+" Remains: "+ StringUtils.niceDuration((vcfs.size()-(i+1))*millisecPerVcf));
					}
				
				}
			w.close();
			w=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(w!=null) try {w.close();} catch(Exception err) {}
			}
		}

	public static void main(final String[] args)	{
		new VcfConcat().instanceMainWithExit(args);
		}

	
	}
