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
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.SAMFileHeader;
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
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFSampleHeaderLine;
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
public class BaseCoverage extends AbstractBaseCov
	{
	private static Logger LOG=Logger.build(BaseCoverage.class).make();

	@Parameter(names={"-Z","--disable-normalization"},description="disable depth normalization _on median")
	protected boolean disable_normalization  = false;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path faidx = null;
	@Parameter(names={"-r","--region","--interval"},description=IntervalParser.OPT_DESC,required=true)
	private String intervalStr=null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private static class SampleInfo {
		String sampleName = null;
		int columnIndex = -1;
		Double medianCoverage = null;
		Double meanCoverage = null;
		final Map<Integer,Integer> coverage2count=new TreeMap<>();
		int count_unmapped=0;
	}
	
	private static class Base {
		int pos;
		int sample_idx;
		int raw_depth;
		double norm_depth;
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
			b.norm_depth = dis.readDouble();
			return b;
			}
		@Override
		public void encode(final DataOutputStream o, final Base b) throws IOException {
			o.writeInt(b.pos);
			o.writeInt(b.sample_idx);
			o.writeInt(b.raw_depth);
			o.writeDouble(b.norm_depth);
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
			
			final List<SampleInfo> samples = new ArrayList<>(bamsIn.size());
			final Map<String, SampleInfo> sample2info = new HashMap<>();
			
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
					final SampleInfo sampleInfo = new SampleInfo();
					final SAMFileHeader header = sr.getFileHeader();
					sampleInfo.sampleName = header.getReadGroups().stream().
							map(S->S.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
					if(sample2info.containsKey(sampleInfo.sampleName)) {
						LOG.error("duplicate sample "+sampleInfo.sampleName+" in bamPath");
						return -1;
						}
					
					sampleInfo.columnIndex = samples.size();
					
					sample2info.put(sampleInfo.sampleName, sampleInfo);
					samples.add(sampleInfo);
				
					final double[] coverage_d = super.getCoverage(sr, queryInterval) ;
					
					sampleInfo.meanCoverage= Arrays.stream(coverage_d).average().orElse(-1);
					sampleInfo.count_unmapped = (int) Arrays.stream(coverage_d).filter(X->X==0).count();
					
					for(final int covx :new int[]{1,5,10,15,20,30,50,100}) {
						sampleInfo.coverage2count.put(covx,  (int) Arrays.stream(coverage_d).filter(X->X>=covx).count());
						}
					
					
					final double[] coverage_norm;
					
					if(!this.disable_normalization) {
						final OptionalDouble od =super.getMedian(coverage_d);
						if(od.isPresent() && od.getAsDouble()>0.0) {
							final double median_cov = od.getAsDouble();
							sampleInfo.medianCoverage = median_cov;
							coverage_norm = super.normalizeOnMedian(coverage_d,median_cov);
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
				for(int array_index=0;array_index< coverage_d.length;++array_index) {
					final Base b = new Base();
            		b.sample_idx = sampleInfo.columnIndex;
            		b.pos = array_index + queryInterval.getStart();
            		b.raw_depth = (int)coverage_d[array_index];
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
			
			try(VariantContextWriter w = this.writingVariantsDelegate.dictionary(dict).open(super.outputFile);
				ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				final VCFInfoHeaderLine infoMeanDP = new VCFInfoHeaderLine("AVG_DP",1,VCFHeaderLineType.Float,"average DP");
				final VCFInfoHeaderLine infoMinDP = new VCFInfoHeaderLine("MIN_DP",1,VCFHeaderLineType.Integer,"min DP");
				final VCFInfoHeaderLine infoMaxDP = new VCFInfoHeaderLine("MAX_DP",1,VCFHeaderLineType.Integer,"max DP");
				final VCFFormatHeaderLine formatMedianDP = new VCFFormatHeaderLine(FORMAT_NORM_DEPH,1,VCFHeaderLineType.Float,"Depth normalized on median");
				final VCFInfoHeaderLine infoMAD = new VCFInfoHeaderLine("MAD",1,VCFHeaderLineType.Float,"median absolute deviation ( https://en.wikipedia.org/wiki/Median_absolute_deviation )");

				for(SampleInfo si: samples) {
					final StringBuilder sb=new StringBuilder("<");
					sb.append("ID=").append(si.sampleName);
					if(si.meanCoverage!=null && si.meanCoverage>=0) sb.append(",MEAN_COV=").append(si.meanCoverage);
					if(si.medianCoverage!=null && si.medianCoverage>=0) sb.append(",MEDIAN_COV=").append(si.medianCoverage);
					sb.append(",UNMAPPED=").append(si.count_unmapped);
					sb.append(",UNMAPPED_PCT=").append((int)(100.0*si.count_unmapped/(double)queryInterval.getLengthOnReference()));
					for(Integer covx: si.coverage2count.keySet()) {
						sb.append(",COV_GE_"+covx+"=").append(si.coverage2count.get(covx));
						sb.append(",COV_GE_"+covx+"_PCT=").append((int)(100.0*si.coverage2count.get(covx)/(double)queryInterval.getLengthOnReference()));
						}
					
					sb.append(">");
					
					final VCFSampleHeaderLine snH = new VCFSampleHeaderLine(sb.toString(), VCFHeaderVersion.VCF4_3);
					metaData.add(snH);
				}
				
				
				metaData.add(infoMeanDP);
				metaData.add(infoMinDP);
				metaData.add(infoMaxDP);
				metaData.add(formatMedianDP);
				metaData.add(infoMAD);
				VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.DEPTH_KEY);
				VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.DEPTH_KEY);
				final VCFHeader header = new VCFHeader(metaData,samples.stream().map(SN->SN.sampleName).collect(Collectors.toList()));
				header.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, header);
				w.writeHeader(header);
				try(CloseableIterator<Base> iter0= sorting.iterator()) {
					final EqualRangeIterator<Base> iter1 = new EqualRangeIterator<>(iter0, (A,B)->A.compare2(B));
					while(iter1.hasNext()) {
						final List<Base> array = iter1.next();
						final Base first = array.get(0);
						final Map<String, Integer> sample2depth = new HashMap<>(samples.size());
						array.stream().forEach(B->sample2depth.put(samples.get(B.sample_idx).sampleName, B.raw_depth));
						final Map<String, Double> sample2medp;
						
						if(!this.disable_normalization) {
							sample2medp = new HashMap<>(samples.size());
							array.stream().filter(B->B.norm_depth>=0.0).forEach(B->sample2medp.put(samples.get(B.sample_idx).sampleName, B.norm_depth));
							}
						else
							{
							sample2medp = Collections.emptyMap();
							}
						
						final List<Genotype> genotypes = new ArrayList<>(samples.size());
						for(final SampleInfo sn:samples) {
							final GenotypeBuilder gb=new GenotypeBuilder(sn.sampleName);
							gb.DP(sample2depth.getOrDefault(sn.sampleName, 0));
							final Double mDP = sample2medp.getOrDefault(sn.sampleName, null);
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
						
						
						vcb.attribute(VCFConstants.DEPTH_KEY, samples.stream().mapToInt(S->sample2depth.getOrDefault(S.sampleName, 0)).sum());
						vcb.attribute(infoMeanDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S.sampleName, 0)).average().orElse(0.0));
						vcb.attribute(infoMinDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S.sampleName, 0)).min().orElse(0));
						vcb.attribute(infoMaxDP.getID(), samples.stream().mapToInt(S->sample2depth.getOrDefault(S.sampleName, 0)).max().orElse(0));
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
