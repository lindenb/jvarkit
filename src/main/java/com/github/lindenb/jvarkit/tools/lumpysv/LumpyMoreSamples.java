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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
BEGIN_DOC

## Example

```
 java -jar dist/lumpyvcf2circos.jar  --minsu 50 -inv -bnb -dup  -o tmp  LumpyExpress.vcf.gz \
  && (cd tmp; /path/to/bin/circos  -outputdir . -conf lumpy.circos.conf  )
```

END_DOC
*/
@Program(name="lumpymoresamples",
		description="TODO",
		keywords={"lumpy","vcf"},
		generate_doc=false
		)
public class LumpyMoreSamples extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpyMoreSamples.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-b","--bam"},description="File containing list of bam, one path per line.",required=true)
	private File bamFileList = null;
	
	
	private int[] confidenceIntervalPos(final VariantContext ctx)
		{
		final List<Integer> L=ctx.getAttributeAsIntList("CIPOS", 0);
		if(L.size()!=2) throw new IllegalArgumentException(L.toString());
		return new int[] {L.get(0).intValue(),L.get(1).intValue()};
		}
	private int[] confidenceIntervalEnd(final VariantContext ctx)
		{
		final List<Integer> L=ctx.getAttributeAsIntList("CIEND", 0);
		if(L.size()!=2) throw new IllegalArgumentException(L.toString());
		return new int[] {L.get(0).intValue(),L.get(1).intValue()};
		}

	private static class SupportingReads
		{		
		int splitReads=0;
		int pairedReads=0;
		}
	
	private List<SAMRecord> extractSupportingReads(
			final VariantContext ctx,
			final String sample,
			final SamReader reader,
			QueryInterval intervals[]
			) {
		intervals = QueryInterval.optimizeIntervals(intervals);
		final List<SAMRecord> L = new ArrayList<>();
		final CloseableIterator<SAMRecord> iter= reader.query(intervals,false);
		while(iter.hasNext())
			{
			final SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getCigar()==null || rec.getCigar().isEmpty()) continue;
			final SAMReadGroupRecord rg = rec.getReadGroup();
			if(rg==null) continue;
			if(!sample.equals(rg.getSample())) continue;
			L.add(rec);
			}
		iter.close();
		return L;
		}
	
	@Override
	public int doWork(final List<String> args) {
	VCFIterator r=null;
	VariantContextWriter vcw=null;
	final Map<String,SamReader> sample2samreaders = new HashMap<>();

	try {
		r = super.openVCFIterator(oneFileOrNull(args));
		
		final VCFHeader headerIn =r.getHeader();
		final SAMSequenceDictionary dict = headerIn.getSequenceDictionary();
		if(dict==null)
			{
			LOG.error(JvarkitException.VcfDictionaryMissing.getMessage("input vcf"));
			return -1;
			}
		
		final SamReaderFactory samReaderFactory = SamReaderFactory.
				makeDefault().
				validationStringency(ValidationStringency.LENIENT);
		IOUtil.slurpLines(this.bamFileList).stream().forEach(F->{
			if(F.trim().isEmpty()) return;
			final SamReader sr = samReaderFactory.open(SamInputResource.of(F));
			final SAMFileHeader samHeader = sr.getFileHeader();
			final SAMSequenceDictionary dict2 = samHeader.getSequenceDictionary();
			if(dict2==null)
				{
				throw new JvarkitException.BamDictionaryMissing(F);
				}
			if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict2)) {
				throw new JvarkitException.DictionariesAreNotTheSame(dict,dict2);
				}
			for(final SAMReadGroupRecord rg: samHeader.getReadGroups())
				{
				final String sample = rg.getSample();
				if(StringUtil.isBlank(sample)) continue;
				final SamReader reader = sample2samreaders.get(sample);
				if(reader==null) {
					sample2samreaders.put(sample,reader);
					}
				else if(reader==sr)
					{
					continue;
					}
				else
					{
					throw new JvarkitException.UserError("more than one sample per bam:"+F);
					}
				}
			});
		final Set<String> inVcfSampleNames = new HashSet<>(headerIn.getSampleNamesInOrder());
		final Set<String> outVcfSampleNames = new HashSet<>(inVcfSampleNames);
		outVcfSampleNames.addAll(sample2samreaders.keySet());
		
		final VCFHeader headerOut =
				new VCFHeader(
						headerIn.getMetaDataInInputOrder(),
						outVcfSampleNames
						);
		
		
		final VCFFormatHeaderLine SU2= new VCFFormatHeaderLine("SU2",1,VCFHeaderLineType.Integer,"Number of pieces of evidence supporting the variant");
		final VCFFormatHeaderLine PE2= new VCFFormatHeaderLine("PE2",1,VCFHeaderLineType.Integer,"Number of split reads supporting the variant");
		final VCFFormatHeaderLine SR2= new VCFFormatHeaderLine("SR2",1,VCFHeaderLineType.Integer,"Number of paired-end reads supporting the variant");
		headerOut.addMetaDataLine(SU2);
		headerOut.addMetaDataLine(PE2);
		headerOut.addMetaDataLine(SR2);
		
		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(headerOut);
		while(r.hasNext())
			{
			final VariantContext ctx = r.next();
			final StructuralVariantType sttype=ctx.getStructuralVariantType();
			if(sttype==null) continue;
			final int tid=dict.getSequenceIndex(ctx.getContig());

			
			final Map<String,Genotype> genotypeMap = new HashMap<>();
			ctx.getGenotypes().stream().forEach(G->genotypeMap.put(G.getSampleName(),G));
			
			for(final String sample: sample2samreaders.keySet())
				{
				final SamReader samReader = sample2samreaders.get(sample);
				final SupportingReads sr = new SupportingReads();
				switch(sttype)
					{
					case DEL:
						{
						int pos = ctx.getStart();
						int ci[]=confidenceIntervalPos(ctx);
						final QueryInterval left = new QueryInterval(tid,pos-ci[0],pos+ci[1]);
						
						int end = ctx.getEnd();
						ci = confidenceIntervalEnd(ctx);
						final QueryInterval right = new QueryInterval(tid,end-ci[0],end+ci[1]);
						
						
						for(final SAMRecord rec: extractSupportingReads(ctx, sample, samReader, new QueryInterval[] {left,right}))
							{
							final Cigar cigar=rec.getCigar();
							if(cigar.isLeftClipped())
								{
								final QueryInterval qi = new QueryInterval(tid,rec.getUnclippedStart(),rec.getStart()-1);
								if(qi.overlaps(left))
									{
									sr.splitReads++;
									if(rec.getReadPairedFlag()) sr.pairedReads++;
									}
								}
							if(cigar.isRightClipped())
								{
								final QueryInterval qi = new QueryInterval(tid,rec.getEnd()+1,rec.getUnclippedEnd());
								if(qi.overlaps(right))
									{
									sr.splitReads++;
									if(rec.getReadPairedFlag()) sr.pairedReads++;
									}
								}
							}
						break;
						}
					default:break;
					}
				final GenotypeBuilder gb;
				if(genotypeMap.containsKey(sample))
					{
					gb=new GenotypeBuilder(genotypeMap.get(sample));
					}
				else
					{
					gb = new GenotypeBuilder(sample,Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
					}
				gb.attribute(SR2.getID(), sr.splitReads);
				gb.attribute(PE2.getID(), sr.pairedReads);
				gb.attribute(SU2.getID(), 0);
				genotypeMap.put(sample, gb.make());
				}
			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			// add missing samples.
			for(final String sampleName:outVcfSampleNames)
				{
				if(genotypeMap.containsKey(sampleName)) continue;
				genotypeMap.put(sampleName, 
						new GenotypeBuilder(sampleName,Arrays.asList(Allele.NO_CALL,Allele.NO_CALL)).
							make());
				}
			vcb.genotypes(genotypeMap.values());
			
			vcw.add(vcb.make());
			}
		r.close();r=null;
		sample2samreaders.values().stream().
			forEach(R->CloserUtil.close(R));
		
		LOG.info("done");
		return 0;
		} 
	catch(final Exception err)
			{
				LOG.error(err);
				return -1;
			}
		finally {
			CloserUtil.close(r);		
		}
	}
	
	
public static void main(String[] args) {
	new LumpyMoreSamples().instanceMainWithExit(args);
}
}
