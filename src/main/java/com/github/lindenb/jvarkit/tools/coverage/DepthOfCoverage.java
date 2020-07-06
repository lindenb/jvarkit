/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.coverage;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.List;
import java.util.Objects;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC

## Example

```
$ java  -jar dist/depthofcoverage.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null  | column -t 

#BAM                       Sample  Contig  Length  Count   Depth
src/test/resources/S1.bam  S1      RF01    3302    25037   7.582374318594791
src/test/resources/S1.bam  S1      RF02    2687    20275   7.545589877186453
src/test/resources/S1.bam  S1      RF03    2592    19583   7.55516975308642
src/test/resources/S1.bam  S1      RF04    2362    17898   7.577476714648603
src/test/resources/S1.bam  S1      RF05    1579    11887   7.528182393920202
src/test/resources/S1.bam  S1      RF06    1356    10201   7.522861356932153
src/test/resources/S1.bam  S1      RF07    1074    8115    7.555865921787709
src/test/resources/S1.bam  S1      RF08    1059    7980    7.5354107648725215
src/test/resources/S1.bam  S1      RF09    1062    7980    7.5141242937853105
src/test/resources/S1.bam  S1      RF10    751     5740    7.6431424766977365
src/test/resources/S1.bam  S1      RF11    666     5037    7.563063063063063
src/test/resources/S1.bam  S1      *       18490   139733  7.557220118983234
src/test/resources/S2.bam  S2      RF01    3302    25030   7.580254391278014
src/test/resources/S2.bam  S2      RF02    2687    20272   7.544473390398213
src/test/resources/S2.bam  S2      RF03    2592    19592   7.558641975308642
src/test/resources/S2.bam  S2      RF04    2362    17916   7.585097375105843
src/test/resources/S2.bam  S2      RF05    1579    11892   7.531348955034832
src/test/resources/S2.bam  S2      RF06    1356    10217   7.534660766961652
src/test/resources/S2.bam  S2      RF07    1074    8112    7.553072625698324

```

END_DOC
 */
@Program(name="depthofcoverage",
	description="A custom 'Depth of Coverage'.",
	keywords={"depth","bam","sam","coverage"},
	creationDate="20190927",
	modificationDate="20200706"
	)
public class DepthOfCoverage extends Launcher
	{
	private static Logger LOG=Logger.build(DepthOfCoverage.class).make();

	
	@Parameter(names={"-R","--reference"},description="For reading CRAM. " + INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"-B","--mask"},description="optional bed containing regions to be MASKED")
	private Path maskBed = null;
	@Parameter(names={"--mapq"},description=" min mapping quality.")
	private int mapping_quality=1;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--auto-mask"},description="Use REFerence sequence to automatically mask bases that are not ATGC")
	private boolean auto_mask = false;
	@Parameter(names={"--skip"},description="Chromosomes to skip (regular expression)")
	private String skipContigExpr = "(NC_007605|hs37d5)";
	@Parameter(names={"--async"},description="use async I/O",hidden=true)
	private boolean asyncIo=false;
	@Parameter(names={"--disable-paired-overlap"},description="Count overlapping bases with mate for paired-end")
	private boolean disable_paired_overlap_flag=false;
	@Parameter(names={"--max-depth"},description="Ignore depth if it is bigger than this value.")
	private int max_depth = 10_000_000;

	
	@Override
	public int doWork(final List<String> args)
		{
		PrintWriter out=null;
		if(this.auto_mask && this.faidx==null) {
			LOG.error("Cannot auto mask if REF is not defined");
			return -1;
			}
		ReferenceSequenceFile referenceSequenceFile=null;
		try
			{
			final Predicate<String> isRejectContig;
			if(!StringUtils.isBlank(this.skipContigExpr)) {
				final Pattern pat = Pattern.compile(this.skipContigExpr);
				isRejectContig = S-> pat.matcher(S).matches();
				}
			else
				{
				isRejectContig = S -> false;
				}	
			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) 
				{
				srf.referenceSequence(this.faidx);
				srf.setUseAsyncIo(this.asyncIo);
				referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
				}
			
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			out.println("#BAM\tSample\tContig\tContig-Length\tMasked-Contig-Length\tCount\tDepth\tMedian");
			
			for(final Path path: IOUtils.unrollPaths(args)) {
				
				try(final SamReader sr = srf.open(path)) {
					if(!sr.hasIndex()) {
						LOG.error("File "+path+" is not indexed.");
						return -1;
					}
					final SAMFileHeader header = sr.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					final Set<String> rejectContigSet = dict.getSequences()
							.stream()
							.map(SSR->SSR.getSequenceName())
							.filter(isRejectContig)
							.collect(Collectors.toSet())
							;
							
					
					
					if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
						LOG.error("file is not sorted on coordinate :"+header.getSortOrder()+" "+path);
						return -1;
						}
					long count_raw_bases = 0L;
					long count_bases = 0L;
					long sum_coverage = 0L;
					final DiscreteMedian<Integer> discreteMedian_wg = new DiscreteMedian<>();
					final String sample = header.getReadGroups().
							stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().orElse(path.toString())
							;
					int coverage[] = null;
					String prevContig = null;
					
					BitSet mask=null;
					final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
					try(CloseableIterator<SAMRecord> iter=sr.iterator()) {
						for(;;)
							{
							final SAMRecord rec = iter.hasNext()?progress.apply(iter.next()):null;
							
							if(rec!=null) {
								if(!SAMRecordDefaultFilter.accept(rec,this.mapping_quality)) continue;
								if(rejectContigSet.contains(rec.getContig())) continue;
								}
							
							if(rec==null || !rec.getContig().equals(prevContig)) {
								if(coverage!=null) {//DUMP
									long count_bases_ctg = 0L;
									long sum_coverage_ctg = 0L;
									final DiscreteMedian<Integer> discreteMedian_ctg = new DiscreteMedian<>();
									
									for(int i=0;i< coverage.length;i++) {
										if(mask.get(i)) continue;
										final int covi = coverage[i];
										if(covi> this.max_depth) continue;
										count_bases_ctg++;
										sum_coverage_ctg += covi;
										discreteMedian_ctg.add(covi);
										}
									out.print(path);
									out.print("\t");
									out.print(sample);
									out.print("\t");
									out.print(prevContig);
									out.print("\t");
									out.print(coverage.length);
									out.print("\t");
									out.print(count_bases_ctg);
									out.print("\t");
									out.print(sum_coverage_ctg);
									out.print("\t");
									if(count_bases_ctg>0) {
										out.printf("%.2f",sum_coverage_ctg/(double)count_bases_ctg);
										}
									else
										{
										out.print("N/A");
										}
									out.print("\t");
									final OptionalDouble median = discreteMedian_ctg.getMedian();
									if(median.isPresent()) {
										out.print(median.getAsDouble());
										}
									else
										{
										out.print("N/A");
										}
									out.println();
									
									count_bases += count_bases_ctg;
									sum_coverage += sum_coverage_ctg;
									count_raw_bases += coverage.length;
									discreteMedian_wg.add(discreteMedian_ctg);
									}
								coverage=null;
								mask=null;
								///
								System.gc();
								if(rec==null) break;
								
								final SAMSequenceRecord ssr = Objects.requireNonNull(dict.getSequence(rec.getContig()));
								coverage = new int[ssr.getSequenceLength()];
								mask = new BitSet(ssr.getSequenceLength());
								if(this.auto_mask && referenceSequenceFile!=null) {
									final byte refSeq[] = Objects.requireNonNull(referenceSequenceFile.getSequence(ssr.getSequenceName())).getBases();
									for(int i=0;i< refSeq.length;i++) {
										if(AcidNucleics.isATGC(refSeq[i])) continue;
										mask.set(i);
									}
								}
								
								/* read mask */
								if(this.maskBed!=null ) {
									final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
									final BedLineCodec codec= new BedLineCodec();
									try(BufferedReader br=IOUtils.openPathForBufferedReading(this.maskBed)) {
										String line;
										while((line=br.readLine())!=null) {
											final BedLine bed = codec.decode(line);
											if(bed==null) continue;
											String ctg = contigNameConverter.apply(bed.getContig());
											if(StringUtils.isBlank(ctg)) continue;
											if(!rec.getContig().equals(ctg)) continue;
											for(int p1=bed.getStart();p1<=bed.getEnd() && p1 <= coverage.length;++p1) {
												mask.set(p1-1);
												}
											}
										}
									}
								
								prevContig=rec.getContig();
								}
							
							int max_end1 = coverage.length;
							
							if(this.disable_paired_overlap_flag && 
								rec.getReadPairedFlag() && 
								!rec.getMateUnmappedFlag() &&
								rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
								rec.getAlignmentStart() < rec.getMateAlignmentStart() &&
								rec.getAlignmentEnd() > rec.getMateAlignmentStart()
								) {
								max_end1 = rec.getMateAlignmentStart() - 1;
								}
							
							for(final AlignmentBlock block:rec.getAlignmentBlocks()) {
								final int pos1=block.getReferenceStart();
								final int len = block.getLength();
								for(int i=0;i< len;i++) {
									if(pos1>0 && pos1 <= max_end1) {
										coverage[pos1-1]++;
										}
									}
								}
							
							}/* end rec */
					
					
						} /* end iter */
					progress.close();
					
					out.print(path);
					out.print("\t");
					out.print(sample);
					out.print("\t");
					out.print(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
					out.print("\t");
					out.print(count_raw_bases);
					out.print("\t");
					out.print(count_bases);
					out.print("\t");
					out.print(sum_coverage);
					out.print("\t");
					if(count_bases>0) {
						out.printf("%.2f",sum_coverage/(double)count_bases);
						}
					else
						{
						out.print("N/A");
						}
					out.print("\t");
					final OptionalDouble median = discreteMedian_wg.getMedian();
					if(median.isPresent()) {
						out.print(median.getAsDouble());
						}
					else
						{
						out.print("N/A");
						}
					out.println();
					}
				}
			out.flush();
			out.close();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(referenceSequenceFile);
			}

		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new DepthOfCoverage().instanceMainWithExit(args);
		}		

}
