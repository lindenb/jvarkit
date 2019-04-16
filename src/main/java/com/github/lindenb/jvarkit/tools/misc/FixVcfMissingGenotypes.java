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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

/**

BEGIN_DOC

## Slow

This tool remains slow because there is a random-access in the bam for each './.' genotype.

You can always try to speed-up things by breaking your VCF in multiple regions and process them in parallel.

## Examples


### Example

```
$ find ~/src/gatk-ui/testdata/ -name "*.bam" > input.list

$ tail -2 input.vcf
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	0/0:0,244,70	0/0:0,199,65	0/0:0,217,68	1/1:69,84,0
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	./.	./.	./.	./.

$ java -jar dist/fixvcfmissinggenotypes.jar -d 50 --fixDP --filtered zz -B input.list input.vcf | tail -2
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=188;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:PL	0/0:48:0,244,70	0/0:63:0,199,65	0/0:53:0,217,68	1/1:24:69,84,0
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=72;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:FT:FXG	./.:48:PASS	0/0:63:zz:1	0/0:53:zz:1	./.:24:PASS
```

### Example


```
$ yourtool-mergingvcf 1.vcf 2.vcf 3.vcf > merged.vcf
$ find ./ -name "*.bam" > bams.list
$  java -jar dist/fixvcfmissinggenotypes.jar -B bams.list < merged.vcf > out.vcf
```

```
$ find DIR1 -name "PREFIX_*final.bam"  | grep -E '(S1|S2|S3|S4)' ) > bams.list

$ find DIR1 -name "PREFIX_*_variations.gatk.annotations.vcf.gz" |\
grep -E '(S1|S2|S3|S4)' |\
xargs perl  vcftools_0.1.12b/perl vcftools_0.1.12b/bin/vcf-merge |\
java -jar dist/fixvcfmissinggenotypes.jar -d 10 -f  bams.list |\
gzip --best > out.vcf.gz

```

## Cited in 

 * "Exome sequencing in genomic regions related to racing performance of Quarter Horses" Pereira, G.L., Malheiros, J.M., Ospina, A.M.T. et al. J Appl Genetics (2019). https://doi.org/10.1007/s13353-019-00483-1
 

### History

 * 2018-11-20 : adding features for structural variants
 * 2017-07-24 : rewrite whole program 
 * 2014: Creation


END_DOC
*/


@Program(name="fixvcfmissinggenotypes",
description="After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then a missing genotype is said hom-ref.",
biostars={119007,263309,276811,302581},
keywords={"sam","bam","vcf","sv","genotype"}
)
public class FixVcfMissingGenotypes extends Launcher
	{
	private static final Logger LOG = Logger.build(FixVcfMissingGenotypes.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--depth"},description="minimal depth before setting a genotype to HOM_REF")
	private int minDepth = 10 ;
	@Parameter(names={"-B","--bams"},description="path of indexed BAM path with read Groups. You can put those paths in a text file having a *.list sufffix")
	private List<String> bamList = new ArrayList<>();
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-T","--tag"},description="FORMAT 'Tag' for fixed genotype")
	private String fixedTag = "FXG";
	@Parameter(names={"--force","-f"},description="[20181120] Update all fields like DP even if the Genotype is called.")
	private boolean forceUpdate=false;
	@Parameter(names={"--filtered"},description="Mark fixed genotypes as FILTERED with this FILTER")
	private String fixedGenotypesAreFiltered=null;
	@Parameter(names={"--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition= SAMRecordPartition.sample;
	
	@ParametersDelegate
	private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

	
	@Override
	protected int doVcfToVcf(final String inputName,final VCFIterator in,final VariantContextWriter out) {
		final List<File> bamFiles=  IOUtils.unrollFiles2018(this.bamList);
		final Map<String,List<SamReader>> sample2bam = new HashMap<>(bamFiles.size());
		final SamReaderFactory srf = super.createSamReaderFactory();
		try {
			final VCFHeader header=in.getHeader();
			for(final File bamFile: bamFiles)
				{
				LOG.info("Reading header for "+bamFile);
				final SamReader reader= srf.open(bamFile);
				if(!reader.hasIndex())
					{
					LOG.error("No BAM index available for "+bamFile);
					return -1;
					}		
				final SAMFileHeader samHeader=reader.getFileHeader();
				for(final SAMReadGroupRecord g:samHeader.getReadGroups())
					{
					final String sample = this.partition.apply(g);
					if(StringUtil.isBlank(sample)) continue;
					if(!header.getSampleNamesInOrder().contains(sample)) continue;
					List<SamReader> readers = sample2bam.get(sample);
					if(readers==null)
						{
						readers=new ArrayList<>();
						sample2bam.put(sample,readers);
						}
					readers.add(reader);
					}
				}
			
			
			final VCFHeader h2=new VCFHeader(header);
			if(h2.getFormatHeaderLine(VCFConstants.DEPTH_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
				}
			if(h2.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
				}
			if(h2.getFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
				}
			final String NDISC_TAG="NDISC";
			final String NPROP_TAG="NOTPROP";
			final String NORPHAN_TAG="NORPHAN";
			final String NSUP_TAG="NSUP";
			final String NCLIPPED="NCLIP";
			if(!StringUtil.isBlank(this.fixedTag)) h2.addMetaDataLine(new VCFFormatHeaderLine(this.fixedTag,1,VCFHeaderLineType.Integer,"Genotype was set as homozygous (min depth ="+this.minDepth+")"));
			h2.addMetaDataLine(new VCFFormatHeaderLine(NDISC_TAG,1,VCFHeaderLineType.Integer,"Number of discordant reads (mate mapped on another contig) in the interval"));
			h2.addMetaDataLine(new VCFFormatHeaderLine(NPROP_TAG,1,VCFHeaderLineType.Integer,"Number of not propertly paired reads in the interval"));
			h2.addMetaDataLine(new VCFFormatHeaderLine(NORPHAN_TAG,1,VCFHeaderLineType.Integer,"Number of orphan reads (mate is not mapped) in the interval"));
			h2.addMetaDataLine(new VCFFormatHeaderLine(NSUP_TAG,1,VCFHeaderLineType.Integer,"Number of reads with supplementary alignments (SA tag)"));
			h2.addMetaDataLine(new VCFFormatHeaderLine(NCLIPPED,1,VCFHeaderLineType.Integer,"Number of clipped reads"));
			
			final short SA_TAG = SAMTag.SA.getBinaryTag();
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().logger(LOG).dictionary(header).build();
			this.recalculator.setHeader(h2);
			out.writeHeader(h2);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.apply(in.next());
				boolean somethingWasChanged=false;
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype genotype = ctx.getGenotype(i);
					final String sample = genotype.getSampleName();
					if(!genotype.isCalled() || this.forceUpdate)
						{
						int depth =0;
						int orphans = 0;
						int notProperlyPaired = 0;
						int dicordantContigs = 0;
						int readWithSupplAligns = 0;
						int clippedReads = 0;
						List<SamReader> samReaders = sample2bam.get(sample);
						if(samReaders==null) samReaders=Collections.emptyList();
						
						for(final SamReader sr: samReaders)
							{
							final SAMRecordIterator iter=sr.query(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false);
							while(iter.hasNext())
								{
								final SAMRecord rec=iter.next();
								if(rec.getReadUnmappedFlag()) continue;
								if(this.filter.filterOut(rec)) continue;
								if(rec.getReadPairedFlag()) {
									if(!rec.getProperPairFlag()) {
										++notProperlyPaired;
										}
									if(rec.getMateUnmappedFlag())
										{
										++orphans;
										}
									else if(!rec.getContig().equals(rec.getMateReferenceName()))
										{
										++dicordantContigs;
										}
									}
								if(rec.getAttribute(SA_TAG)!=null) readWithSupplAligns++;
								final SAMReadGroupRecord rg=rec.getReadGroup();
								if(!sample.equals(this.partition.apply(rg))) continue;
								final Cigar cigar=rec.getCigar();
								if(cigar==null || cigar.isEmpty()) continue;
								if(cigar.isClipped())
									{
									clippedReads ++;
									}
								int refPos = rec.getAlignmentStart();
								for(final CigarElement ce:cigar.getCigarElements())
									{
									if( refPos > ctx.getEnd() ) break;
									if(!ce.getOperator().consumesReferenceBases()) continue;
									if( ce.getOperator().consumesReadBases())
										{
										for(int n=0;n< ce.getLength();++n )
											{
											if( refPos+n < ctx.getStart() ) continue;
											if( refPos+n > ctx.getEnd()) break;
											depth++;
											}
										}
									refPos+= ce.getLength();
									}
								}
							iter.close();
							}
						depth/= 1+ctx.getEnd()-ctx.getStart();
						
						final GenotypeBuilder  gb = new GenotypeBuilder(genotype);
						somethingWasChanged=true;
						gb.DP(depth);
						if(ctx.getStart()<ctx.getEnd()) {
							gb.attribute(NDISC_TAG, dicordantContigs);
							gb.attribute(NPROP_TAG, notProperlyPaired);
							gb.attribute(NORPHAN_TAG, orphans);
							gb.attribute(NSUP_TAG, readWithSupplAligns);
							gb.attribute(NCLIPPED, clippedReads);
							}
						if(!genotype.isCalled() && depth>= this.minDepth) {
							gb.alleles(Arrays.asList(ctx.getReference(),ctx.getReference()));
							if(!StringUtil.isBlank(this.fixedTag)) gb.attribute(this.fixedTag, 1);
							if(!StringUtil.isBlank(this.fixedGenotypesAreFiltered)) gb.filter(this.fixedGenotypesAreFiltered);
							}
						genotype = gb.make();
						}
					genotypes.add(genotype);
					}//end of for-each genotype
				
				if(somethingWasChanged)
					{
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.genotypes(genotypes);
					out.add(this.recalculator.apply(vcb.make()));
					}
				else
					{
					out.add(ctx);
					}
				}
			progress.close();

			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{	
			sample2bam.values().stream().flatMap(L->L.stream()).forEach(R->CloserUtil.close(R));
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	public static void main(final String[] args) {
		new FixVcfMissingGenotypes().instanceMainWithExit(args);
	}

}
