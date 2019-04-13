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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**

BEGIN_DOC

## input

input is a set of bam files or one file with the suffix '.list' containing the path to the bams.

## Example

```
$ java -jar dist/samtranslocations.jar src/test/resources/HG02260.transloc.chr9.14.bam | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG02260
9	137230996	9:137230996:14:79839048	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=14;DP=18;POS2=79839048;STDDEV_POS1=120;STDDEV_POS2=187;SVTYPE=BND	GT:DP:SR	0/1:18:5,13,13,5
14	79839131	14:79839131:9:137230969	N	<TRANSLOC>	17	.	AC=1;AF=0.500;AN=2;CHROM2=9;DP=17;POS2=137230969;STDDEV_POS1=153;STDDEV_POS2=153;SVTYPE=BND	GT:DP:SR	0/1:17:12,5,5,12
```

## History

* 2018-09-18 :  rewriting
* 2017-12-13 :  refactoring for balanced translocation.

END_DOC
*/
@Program(name="samtranslocations",
	description="Explore balanced translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam","sv","translocation"},
	modificationDate="20190329"
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.build(SamTranslocations.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-d","--distance"},
			description="Max distance between two read to test if they both end at the same ~ position. " +DistanceParser.OPT_DESCRIPTION,
			converter=DistanceParser.StringConverter.class,
			splitter=NoSplitter.class
			)
	private int fuzzy_distance = 1_000;
	@Parameter(names={"-m","--min"},description="Min number of events to validate the translocation")
	private int min_number_of_events=3;
	@Parameter(names={"-B","--bed"},description="Optional BED file. SV should overlap this bed.")
	private Path bedFile = null;
	@Parameter(names={"--chrom-regex"},description="Only consider the chromosomes matching the following regular expression.")
	private String chromRegex="(chr)?([1-9][0-9]*|X|Y)";
	@Parameter(names={"--low-complexity"},description="zone of low complexity (=many discordant reads). In those region the software will 'freeze'. "
			+ "If the number of buffered reads * this number per sample then clear the buffer")
	private long max_complexity_buffer = 2_000L;
	@Parameter(names={"--max-sa"},description="ignore reads having more that SA:Z supplementary alignments.")
	private long max_sa_attribute = 4;
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 0;

	/* for the data I tested there was a bug in the sorting (reads not sorted on FLAG see htsjdk.samtools.SAMRecordCoordinateComparator)
	 but we just need fileOrderCompare */
	
	
	@SuppressWarnings("serial")
	private final Comparator<SAMRecord> coordinateComparator = new htsjdk.samtools.SAMRecordCoordinateComparator() {
		@Override
		public int compare(final SAMRecord samRecord1,final SAMRecord samRecord2)
			{
			return super.fileOrderCompare(samRecord1, samRecord2);
			}
		};
	
	private CloseableIterator<SAMRecord> makeIterator(final SamReader sr) {
			try {
			CloseableIterator<SAMRecord> iter1;
			final SAMFileHeader header=sr.getFileHeader();
			if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				throw new JvarkitException.BamBadSortOrder(SAMFileHeader.SortOrder.coordinate, header.getSortOrder());
				}
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final Pattern pat = StringUtils.isBlank(this.chromRegex)?null:Pattern.compile(this.chromRegex);
	
			
			final boolean allowedContigs[] = new boolean[dict.size()];
			Arrays.fill(allowedContigs, true);
			if(pat!=null) {
				dict.getSequences().stream().
					filter(SR->!pat.matcher(SR.getSequenceName()).matches()).
					forEach(SR->allowedContigs[SR.getSequenceIndex()]=false);
				}

			
			if(this.bedFile==null) {
				iter1 = sr.iterator();
				}
			else
				{
				if(!sr.hasIndex()) {
					throw new RuntimeIOException("SamReader "+sr.getResourceDescription()+" is not indexed");
				}
				final BedLineCodec bedCodec=new BedLineCodec();
				final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
				final List<QueryInterval> queryIntervals = new ArrayList<>();
				try(BufferedReader br=com.github.lindenb.jvarkit.io.IOUtils.openPathForBufferedReading(this.bedFile)) {
					br.lines().
					filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
					map(line->bedCodec.decode(line)).
					filter(B->B!=null).
					map(B->B.toInterval()).
					filter(L->L.getStart()<L.getEnd()).
					forEach(B->{
						final String c = contigNameConverter.apply(B.getContig());
						if(StringUtils.isBlank(c)) return;
						final int tid = dict.getSequenceIndex(c);
						if(tid<0) return;
						if(!allowedContigs[tid]) return;
						final QueryInterval qi = new QueryInterval(tid,B.getStart(),B.getEnd());
						queryIntervals.add(qi);							
						});	
					}
				final QueryInterval array[] = QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
				iter1=sr.queryOverlapping(array);
				}
			
	
			
			final short SA = SAMTag.SA.getBinaryTag();
			
			final FilterIterator<SAMRecord> fsi=new FilterIterator<>(iter1, 
					(SR)->{
						if(SR.getReadUnmappedFlag())  return false;
						if(SR.getDuplicateReadFlag()) return false;
						if(SR.isSecondaryOrSupplementary()) return false;
						if(SR.getMappingQuality()<=min_mapq) return false;
						final SAMReadGroupRecord rg=SR.getReadGroup();
						if(rg==null) return false;
						final String sn  = rg.getSample();
						if(StringUtils.isBlank(sn)) return false;
						
						if(max_sa_attribute>0 &&
							SR.getAttribute(SA)!=null &&
							SAMUtils.getOtherCanonicalAlignments(SR).
							stream().
							filter(R->allowedContigs[R.getReferenceIndex()]).
							count() > max_sa_attribute) {
							return false;
							}
						
						if(SR.getReadPairedFlag() &&
							!SR.getMateUnmappedFlag() &&
							allowedContigs[SR.getReferenceIndex()] &&
							!SR.getReferenceIndex().equals(SR.getMateReferenceIndex()) &&
							allowedContigs[SR.getMateReferenceIndex()]
							) {
							return true;
							}
						
						return false;
						}
					);
			return fsi;
			} 
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	private void putbackInBuffer(final List<SAMRecord> candidates,final List<SAMRecord> buffer) {
		buffer.addAll(candidates.subList(1, candidates.size()));//remove first element
		Collections.sort(buffer,this.coordinateComparator);
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.fuzzy_distance<=0) {
			LOG.error("fuzzy_distance <=0 ("+fuzzy_distance+")");
			return -1;
		}
		final List<SamReader> samReaders = new ArrayList<>();
		final List<CloseableIterator<SAMRecord>> samRecordIterators = new ArrayList<>();
		VariantContextWriter out = null;
		try {
			final List<Path> inputBams = IOUtils.unrollPaths(args);
			final SamReaderFactory srf = super.createSamReaderFactory();
		
			if(inputBams.isEmpty()) {
				LOG.error("No bam was defined");
				return -1;
				}
			
			samReaders.addAll( inputBams.stream().map(P->srf.open(P)).collect(Collectors.toList()));
			final SAMFileHeader header = new SamFileHeaderMerger(SortOrder.coordinate,
					samReaders.stream().map(SR->SR.getFileHeader()).collect(Collectors.toList()),
					false).getMergedHeader();
			samRecordIterators.addAll(samReaders.stream().map(S->makeIterator(S)).collect(Collectors.toList()));
			final CloseableIterator<SAMRecord>	iter = new MergingIterator<>(
					this.coordinateComparator,
					samRecordIterators
					);
				
				
			final Set<String> sampleNames = header.getReadGroups().
				stream().
				map(RG->RG.getSample()).
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toSet());	
			
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(header);
			
			final boolean allowedContigs[] = new boolean[refDict.size()];
			Arrays.fill(allowedContigs, true);
			if(!StringUtils.isBlank(this.chromRegex)) {
				final Pattern pat = Pattern.compile(this.chromRegex);
				refDict.getSequences().stream().
					filter(SR->!pat.matcher(SR.getSequenceName()).matches()).
					forEach(SR->allowedContigs[SR.getSequenceIndex()]=false);
				}
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
					newInstance().
					dictionary(refDict).
					logger(LOG).
					build();
			
			
			
			
			final Set<VCFHeaderLine> metaData=new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
			metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));
			final VCFFormatHeaderLine supportingReadsFormat = new VCFFormatHeaderLine("SR",
					4,
					VCFHeaderLineType.Integer,
					"Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse"
					);
			metaData.add(supportingReadsFormat);
			final VCFInfoHeaderLine stdDevContig1Info = new VCFInfoHeaderLine("STDDEV_POS1",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 1"
					);
			metaData.add(stdDevContig1Info);
			final VCFInfoHeaderLine stdDevContig2Info = new VCFInfoHeaderLine("STDDEV_POS2",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 2"
					);
			metaData.add(stdDevContig2Info);
			final VCFInfoHeaderLine chrom2Info = new VCFInfoHeaderLine("CHROM2",
					1,
					VCFHeaderLineType.String,
					"other chromosome"
					);
			metaData.add(chrom2Info);
			final VCFInfoHeaderLine pos2Info = new VCFInfoHeaderLine("POS2",
					1,
					VCFHeaderLineType.Integer,
					"other position"
					);
			metaData.add(pos2Info);
			final VCFInfoHeaderLine highQualSampleInfo = new VCFInfoHeaderLine("highQualSamples",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"High Quality Samples"
					);
			metaData.add(highQualSampleInfo);
			final VCFInfoHeaderLine lowQualSampleInfo = new VCFInfoHeaderLine("lowQualSamples",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Low Quality Samples"
					);
			metaData.add(lowQualSampleInfo);
			final VCFInfoHeaderLine maxEventInfo = new VCFInfoHeaderLine("MAX_DP",
					1,
					VCFHeaderLineType.Integer,
					"Max number of event per sample"
					);
			metaData.add(maxEventInfo);


			final VCFHeader vcfHeader= new VCFHeader(metaData, sampleNames);
			vcfHeader.setSequenceDictionary(refDict);
			JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
			
			out = VCFUtils.createVariantContextWriterToPath(this.outputFile);
			out.writeHeader(vcfHeader);
			
			
			final Function<SAMRecord, Integer> mateEnd = REC->SAMUtils.getMateCigar(REC)!=null?
					SAMUtils.getMateAlignmentEnd(REC):
					REC.getMateAlignmentStart()
					;
					
			final Allele REF=Allele.create("N", true);
			final Allele ALT=Allele.create("<TRANSLOC>", false);
			/** buffer or reads */
			final LinkedList<SAMRecord> buffer= new LinkedList<>();

			int prev_tid=-1;
			for(;;)
				{
				final SAMRecord rec;
				if(!buffer.isEmpty()) {
					rec= buffer.pollFirst();
					}
				else if(iter.hasNext()) {
					rec= progress.apply(iter.next());
					}
				else
					{
					rec=null;
					}
			
				if(rec==null || (rec.getReferenceIndex()!=prev_tid))
					{
					final int final_prev_tid = prev_tid;
					buffer.removeIf(SR->SR.getReferenceIndex()<=final_prev_tid);
					if(rec==null) break;
					prev_tid= rec.getReferenceIndex();
					}
				buffer.removeIf(SR->SR.getUnclippedEnd()<rec.getUnclippedStart());

				if(rec.getReadNegativeStrandFlag()) continue;
				
				final List<SAMRecord> candidates = new ArrayList<>();
								
				candidates.add(rec);
				
				if(rec.getReadPairedFlag() && 
					!rec.getReadNegativeStrandFlag() &&
					!rec.getMateUnmappedFlag() &&
					!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())
					)
					{
					long end = rec.getUnclippedEnd() + this.fuzzy_distance;
					
					int buffer_index = 0;
					while(buffer_index < buffer.size())
						{
						final SAMRecord rec2 = buffer.get(buffer_index);
						
						if(!rec2.getReferenceIndex().equals(rec.getReferenceIndex())) {
							break;
							}
						if(rec2.getAlignmentStart() > end) {//not unclipped to avoid side effect
							break;
							}
						if(rec2.getMateUnmappedFlag()) {
							buffer_index++;
							continue;
							}
						if(!rec2.getMateReferenceIndex().equals(rec.getMateReferenceIndex())) {
							buffer_index++;
							continue;
							}
						final int mate1 = rec.getMateNegativeStrandFlag()?
									rec.getMateAlignmentStart():
									mateEnd.apply(rec)
									;
						final int mate2 = rec2.getMateNegativeStrandFlag()?
								rec2.getMateAlignmentStart():
								mateEnd.apply(rec2)
								;

						if(Math.abs(mate1-mate2) > this.fuzzy_distance) {
							buffer_index++;
							continue;
							}
						buffer.remove(buffer_index);
						candidates.add(rec2);						
						}
					
					
					while(iter.hasNext())
						{
						final SAMRecord rec2 = iter.next();
						if(rec2==null) break;

						if(!rec2.getReferenceIndex().equals(rec.getReferenceIndex())) {
							buffer.add(rec2);
							break;
							}
			
						if(rec2.getAlignmentStart() > end) {//not unclipped to avoid side effect
							buffer.add(rec2);
							break;
							}
						if(rec2.getMateUnmappedFlag()) {
							buffer.add(rec2);
							continue;
							}
						if(!rec2.getMateReferenceIndex().equals(rec.getMateReferenceIndex())) {
							buffer.add(rec2);
							continue;
							}
						final int mate1 = rec.getMateNegativeStrandFlag()?
									rec.getMateAlignmentStart():
									mateEnd.apply(rec)
									;
						final int mate2 = rec2.getMateNegativeStrandFlag()?
								rec2.getMateAlignmentStart():
								mateEnd.apply(rec2)
								;

						if(Math.abs(mate1-mate2) > this.fuzzy_distance) {
							buffer.add(rec2);
							continue;
							}
						candidates.add(rec2);
						}
					}
				else
					{
					continue;
					}
				if((long)buffer.size()> this.max_complexity_buffer*samReaders.size())
					{
					//final long end = rec.getUnclippedEnd();
					//buffer.remove(0);
					//buffer.removeIf(R->R.getReferenceIndex().equals(rec.getReferenceIndex()) && R.getStart()<=end);
					buffer.clear();
					LOG.warn("zone of low complexity: clearing buffer near "+ rec.getContig()+":"+rec.getStart());
					continue;
					}
				if(candidates.isEmpty() || candidates.size()<this.min_number_of_events) continue;
				

				//check in both sides
				final int count_plus =  (int)candidates.stream().filter(SR->!SR.getReadNegativeStrandFlag()).count();
				final int count_minus = (int)candidates.stream().filter(SR-> SR.getReadNegativeStrandFlag()).count();
				if(count_plus<this.min_number_of_events || count_minus<this.min_number_of_events) {
					putbackInBuffer(candidates, buffer);
					continue;
				}
				
				
				final int pos_ctg1=(int)candidates.stream().mapToInt(
						SR->SR.getReadNegativeStrandFlag()?SR.getAlignmentStart():SR.getAlignmentEnd()).
						average().orElse(-1);
				if(pos_ctg1<1) {
					putbackInBuffer(candidates, buffer);
					continue;
				}
				final int stddev_ctg1 = (int)candidates.stream().mapToInt(
						SR->SR.getReadNegativeStrandFlag()?SR.getAlignmentStart():SR.getAlignmentEnd()).
						map(X->Math.abs(X-pos_ctg1)).
						average().
						orElse(0.0);
				
				final int pos_ctg2=(int)candidates.stream().
						mapToInt(SR->SR.getMateNegativeStrandFlag()?SR.getMateAlignmentStart():mateEnd.apply(SR)).
						average().orElse(-1);
				if(pos_ctg2<1) {
					putbackInBuffer(candidates, buffer);
					continue;
					}
				final int stddev_ctg2 = (int)candidates.stream().
						mapToInt(SR->SR.getMateNegativeStrandFlag()?SR.getMateAlignmentStart():mateEnd.apply(SR)).
						map(X->Math.abs(X-pos_ctg2)).
						average().
						orElse(0.0);
				
				final List<String> lowQualSamples=new ArrayList<>();
				final List<String> hiQualSamples=new ArrayList<>();
				final VariantContextBuilder vcb=new VariantContextBuilder(null,
						rec.getContig(), 
						pos_ctg1,
						pos_ctg1,
						Arrays.asList(REF,ALT)
						);
				final List<Genotype> genotypes = new ArrayList<>(sampleNames.size());
				int max_dp = 0;
				for(final String sample:sampleNames) {
					final List<SAMRecord> readSample = candidates.stream().
							filter(SR->sample.equals(SR.getReadGroup().getSample())).collect(Collectors.toList());
					if(readSample.isEmpty())
						{
						genotypes.add(GenotypeBuilder.createMissing(sample, 2));
						}
					else
						{
						final GenotypeBuilder gb=new GenotypeBuilder(sample,Arrays.asList(REF,ALT));
						final int sn_contig1_count_plus =  (int)readSample.stream().filter(SR->!SR.getReadNegativeStrandFlag()).count();
						final int sn_contig1_count_minus = (int)readSample.stream().filter(SR-> SR.getReadNegativeStrandFlag()).count();
						final int sn_contig2_count_plus =  (int)readSample.stream().filter(SR->!SR.getMateNegativeStrandFlag()).count();
						final int sn_contig2_count_minus = (int)readSample.stream().filter(SR-> SR.getMateNegativeStrandFlag()).count();
						
						max_dp=Math.max(max_dp, readSample.size());
						
						gb.DP(readSample.size());
						gb.attribute(supportingReadsFormat.getID(),
							new int[] {
								sn_contig1_count_plus,
								sn_contig1_count_minus,
								sn_contig2_count_plus,
								sn_contig2_count_minus
								}
							);
						
						genotypes.add(gb.make());

						if(sn_contig1_count_plus>this.min_number_of_events && sn_contig1_count_minus>this.min_number_of_events &&
						   sn_contig2_count_plus>this.min_number_of_events && sn_contig2_count_minus>this.min_number_of_events)  {
							hiQualSamples.add(sample);
							}
						else
							{
							lowQualSamples.add(sample);
							}

						}
					
					}
				if(hiQualSamples.isEmpty()) {
					putbackInBuffer(candidates, buffer);
					continue;
				}
				vcb.id(rec.getReferenceName()+":"+pos_ctg1+":"+rec.getMateReferenceName()+":"+pos_ctg2);
				vcb.attribute(stdDevContig1Info.getID(), stddev_ctg1);
				vcb.attribute(stdDevContig2Info.getID(), stddev_ctg2);
				vcb.attribute(chrom2Info.getID(), rec.getMateReferenceName());
				vcb.attribute(pos2Info.getID(),pos_ctg2);
				vcb.attribute(maxEventInfo.getID(),max_dp);
				
				
				
				vcb.attribute(VCFConstants.SVTYPE, StructuralVariantType.BND.name());
				vcb.attribute(VCFConstants.DEPTH_KEY,candidates.size());
				vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,(hiQualSamples.size()+lowQualSamples.size()));
				vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,sampleNames.size()*2);
				vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,(hiQualSamples.size()+lowQualSamples.size())/(sampleNames.size()*2.0));
				if(!lowQualSamples.isEmpty()) vcb.attribute(lowQualSampleInfo.getID(), lowQualSamples);
				if(!hiQualSamples.isEmpty()) vcb.attribute(highQualSampleInfo.getID(), hiQualSamples);
				
				vcb.genotypes(genotypes);
				vcb.log10PError(candidates.size()/-10.0);
				vcb.alleles(Arrays.asList(REF,ALT));
				
				out.add(vcb.make());
				
				
				}
			progress.close();
			iter.close();
			samRecordIterators.forEach(S->S.close());
			samRecordIterators.clear();
			samReaders.forEach(CloserUtil::close);
			samReaders.clear();
			out.close();
			out=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			samRecordIterators.forEach(S->S.close());
			samReaders.forEach(S->CloserUtil.close(S));
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
