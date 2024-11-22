/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.basecoverage;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StopWatch;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Input

input is a set of bam files or one file with the '.list' suffix containing the path to the bams

## Example:

```
$ java -jar dist/basecoverage.jar --region "RF03:491-500" -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 

lindenb@asimov:~/src/jvarkit$ java -jar dist/basecoverage.jar --bed jeter.bed -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 

##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##basecoverage.meta=compilation:20220420110509 githash:afbe74ab2 htsjdk:2.24.1 date:20220420110554 cmd:--bed jeter.bed -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
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
RF03	491	.	N	.	.	.	DP=28	DP	3	6	6	8	5
RF03	492	.	N	.	.	.	DP=27	DP	3	5	5	8	6
RF03	493	.	N	.	.	.	DP=27	DP	3	5	5	8	6
RF03	494	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	495	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	496	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	497	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	498	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	499	.	N	.	.	.	DP=25	DP	3	4	4	9	5
RF03	500	.	N	.	.	.	DP=23	DP	3	3	3	9	5


```

END_DOC
 */
@Program(name="basecoverage",
	description="'Depth of Coverage' per base.",
	keywords={"depth","bam","sam","coverage","vcf"},
	creationDate="20220420",
	modificationDate="20241122",
	jvarkit_amalgamion =  true,
	menu="BAM Manipulation"
	)
public class BaseCoverage extends Launcher
	{
	private static Logger LOG=Logger.build(BaseCoverage.class).make();

	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path faidx = null;
	@Parameter(names={"-r","--region","--interval"},description=IntervalParser.OPT_DESC,required=true)
	private String intervalStr=null;
	@Parameter(names={"--mapq"},description=" min mapping quality.")
	private int mapping_quality=1;
	@Parameter(names={"-Z","--disable-normalization"},description="disable depth normalization _on median")
	private boolean disable_normalization  = false;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	private static class Base {
		int pos;
		int sample_idx;
		int raw_depth;
		float norm_depth;
		int compare2(final Base o) {
			return Integer.compare(this.pos, o.pos);
		}	
		int compare1(final Base o) {
			int i = compare2(o);
			if(i!=0) return i;
			i = Integer.compare(this.sample_idx, o.sample_idx);
			return i;
		}
	}
	
	private static class BaseCodec extends AbstractDataCodec<Base> {
		@Override
		public Base decode(DataInputStream dis) throws IOException {
			final Base b = new Base();
			try {
				b.pos =dis.readInt();
				}
			catch(final EOFException err) {
				return null;
				}
			b.sample_idx = dis.readInt();
			b.raw_depth = dis.readInt();
			b.norm_depth = dis.readFloat();
			return b;
			}
		@Override
		public void encode(final DataOutputStream o, final Base b) throws IOException {
			o.writeInt(b.pos);
			o.writeInt(b.sample_idx);
			o.writeInt(b.raw_depth);
			o.writeFloat(b.norm_depth);
			}
		@Override
		public BaseCodec clone() {
			return new BaseCodec();
			}
		}

	
	
	@Override
	public int doWork(final List<String> args)
		{
		SortingCollection<Base>	sorting=null;
		try
			{
			final List<Path> bamsIn = IOUtils.unrollPaths(args);
			if(bamsIn.isEmpty()) {
				LOG.error("input is empty.");
				return -1;
				}
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);

			final Locatable queryInterval  =	new IntervalParser(dict).
					enableWholeContig().
					apply(this.intervalStr).
					orElseThrow(()->new IllegalArgumentException("cannot parse"+this.intervalStr));
			
			
			final SamReaderFactory srf = super.createSamReaderFactory().
					referenceSequence(this.faidx);
			
			final List<String> samples = new ArrayList<>(bamsIn.size());
			final Map<String, Integer> sample2idx = new HashMap<>();
			
			sorting = SortingCollection.newInstance(
					Base.class,
					new BaseCodec(),
					(A,B)->A.compare1(B),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorting.setDestructiveIteration(true);
			
			long millisec_per_samples= 0L;
			int sam_idx=0;
			for(final Path bamPath: bamsIn) {
				final StopWatch stopWatch = new StopWatch();
				stopWatch.start();
				IOUtil.assertFileIsReadable(bamPath);
				LOG.info("scanning "+bamPath+" "+(++sam_idx)+"/" + bamsIn.size());
				try(SamReader sr = srf.open(bamPath)) {
					if(!sr.hasIndex()) {
						throw new IllegalArgumentException("bam "+bamPath+" is not indexed");
						}
					final SAMFileHeader header = sr.getFileHeader();
					final String sn = header.getReadGroups().stream().
							map(S->S.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
					if(sample2idx.containsKey(sn)) {
						LOG.error("duplicate sample "+sn+" in bamPath");
						return -1;
						}
					final int sample_idx = samples.size();
					sample2idx.put(sn, sample_idx);
					samples.add(sn);
					System.gc();
					final int[] coverage_int = new int[queryInterval.getLengthOnReference()];
					Arrays.fill(coverage_int,0);
		

					try(CloseableIterator<SAMRecord> it= sr.queryOverlapping(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd())) {
						while(it.hasNext()) {
							final SAMRecord rec =it.next();
							if(!SAMRecordDefaultFilter.accept(rec, this.mapping_quality)) continue;
							for(AlignmentBlock ab : rec.getAlignmentBlocks()) {
								for(int x=0;x < ab.getLength();++x) {
									final int array_index = (ab.getReferenceStart()+x) - queryInterval.getStart();
									if(array_index<0 ) continue;
									if(array_index>=coverage_int.length) break;
									coverage_int[array_index]++;
									}
								}
							}
						}
					
					final float[] coverage_norm;
					
					if(!this.disable_normalization) {
						final DiscreteMedian<Integer> dm = new DiscreteMedian<>();
						for(int array_index=0;array_index< coverage_int.length;++array_index) {
							dm.add(coverage_int[array_index]);
							}
						final OptionalDouble od = dm.getMedian();
						if(od.isPresent() && od.getAsDouble()>0.0) {
							final double median_cov = od.getAsDouble();
							coverage_norm = new float[coverage_int.length];	
							for(int array_index=0;array_index< coverage_int.length;++array_index) {
								coverage_norm[array_index]=(float)(coverage_int[array_index]/median_cov);
								}
							}
						else
							{
							coverage_norm = null;
							}
						}
					else
						{
						coverage_norm  =null;
						}
				for(int array_index=0;array_index< coverage_int.length;++array_index) {
					final Base b = new Base();
            		b.sample_idx = sample_idx;
            		b.pos = array_index + queryInterval.getStart();
            		b.raw_depth = coverage_int[array_index];
            		b.norm_depth = (coverage_norm==null?-1f:coverage_norm[array_index]);
            		sorting.add(b);
					}						
				// end SAMReader
				stopWatch.stop();
				millisec_per_samples+= stopWatch.getElapsedTime();
				LOG.info("That took "+StringUtils.niceDuration(stopWatch.getElapsedTime())+
						". Elapsed:"+
						StringUtils.niceDuration(millisec_per_samples)+" Remains:"+
						StringUtils.niceDuration((long)((millisec_per_samples/(double)sam_idx))*(bamsIn.size()-sam_idx)));
				}//end samReader
			} // end for each bam
			sorting.doneAdding();
			
			try(VariantContextWriter w = this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
				ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				final VCFInfoHeaderLine infoMeanDP = new VCFInfoHeaderLine("AVG_DP",1,VCFHeaderLineType.Float,"average DP");
				final VCFInfoHeaderLine infoMinDP = new VCFInfoHeaderLine("MIN_DP",1,VCFHeaderLineType.Integer,"min DP");
				final VCFInfoHeaderLine infoMaxDP = new VCFInfoHeaderLine("MAX_DP",1,VCFHeaderLineType.Integer,"max DP");
				final VCFFormatHeaderLine formatMedianDP = new VCFFormatHeaderLine("MD",1,VCFHeaderLineType.Float,"Depth normalized on median");
				final VCFInfoHeaderLine infoMAD = new VCFInfoHeaderLine("MAD",1,VCFHeaderLineType.Float,"median absolute deviation ( https://en.wikipedia.org/wiki/Median_absolute_deviation )");

				metaData.add(infoMeanDP);
				metaData.add(infoMinDP);
				metaData.add(infoMaxDP);
				metaData.add(formatMedianDP);
				metaData.add(infoMAD);
				VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.DEPTH_KEY);
				VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.DEPTH_KEY);
				final VCFHeader header = new VCFHeader(metaData,samples);
				header.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, header);
				w.writeHeader(header);
				try(CloseableIterator<Base> iter0= sorting.iterator()) {
					final EqualRangeIterator<Base> iter1 = new EqualRangeIterator<>(iter0, (A,B)->A.compare2(B));
					while(iter1.hasNext()) {
						final List<Base> array = iter1.next();
						final Base first = array.get(0);
						final Map<String, Integer> sample2depth = new HashMap<>(samples.size());
						array.stream().forEach(B->sample2depth.put(samples.get(B.sample_idx), B.raw_depth));
						final Map<String, Float> sample2medp;
						
						if(!this.disable_normalization) {
							sample2medp = new HashMap<>(samples.size());
							array.stream().filter(B->B.norm_depth>=0f).forEach(B->sample2medp.put(samples.get(B.sample_idx), B.norm_depth));
							}
						else
							{
							sample2medp = Collections.emptyMap();
							}
						
						final List<Genotype> genotypes = new ArrayList<>(samples.size());
						for(final String sn:samples) {
							final GenotypeBuilder gb=new GenotypeBuilder(sn);
							gb.DP(sample2depth.getOrDefault(sn, 0));
							final Float mDP = sample2medp.getOrDefault(sn, null);
							if(mDP!=null) gb.attribute(formatMedianDP.getID(),mDP);
							genotypes.add(gb.make());
							}
						final Allele ref_allele = Allele.create(fasta.getSubsequenceAt(queryInterval.getContig(), first.pos, first.pos).getBases(), true);
						final VariantContextBuilder vcb = new VariantContextBuilder(
								null,
								queryInterval.getContig(),
								first.pos,
								first.pos,
								Collections.singletonList(ref_allele)
								);
						vcb.genotypes(genotypes);
						
						if(!this.disable_normalization) {
							final double[] meds = array.stream().mapToDouble(B->B.norm_depth).filter(V->V>=0f).sorted().toArray();
							if(meds.length>0) {
								final double median = new Median().evaluate(meds);
								for(int x=0;x < meds.length;++x) {
									meds[x] = median - meds[x];
									}
								Arrays.sort(meds);
								final double mad_value = new Median().evaluate(meds);
								vcb.attribute(infoMAD.getID(), mad_value);
								}
							}
						
						
						vcb.attribute(VCFConstants.DEPTH_KEY, samples.stream().mapToInt(S->sample2depth.getOrDefault(S, 0)).sum());
						vcb.attribute(infoMeanDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S, 0)).average().orElse(0.0));
						vcb.attribute(infoMinDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S, 0)).min().orElse(0));
						vcb.attribute(infoMaxDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S, 0)).max().orElse(0));
						w.add(vcb.make());
						}//end iter0
					iter1.close();
					} //end sorting iter
				}//end variant context writer
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}

	public static void main(final String[] args)
		{
		new BaseCoverage().instanceMainWithExit(args);
		}		

}
