package com.github.lindenb.jvarkit.tools.coverage;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
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
	description="Depth of Coverage",
	keywords={"depth","bam","sam","coverage"}
	)
public class DepthOfCoverage extends Launcher
	{
	private static Logger LOG=Logger.build(DepthOfCoverage.class).make();

	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-B","--mask"},description="bed containing regions to be MASKED")
	private Path maskBed = null;
	@Parameter(names={"--mapq"},description="mapping quality.")
	private int mapping_quality=0;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;

	@Override
	public int doWork(List<String> args)
		{
		SamReader sfr=null;
		PrintWriter out=null;
		try
			{
			final SamReaderFactory srf = super.createSamReaderFactory();
			srf.referenceSequence(this.faidx);
			
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			out.println("#BAM\tSample\tContig\tLength\tCount\tDepth");
			
			for(final Path path: IOUtils.unrollPaths(args)) {
				
				try(final SamReader sr = srf.open(path)) {
					final SAMFileHeader header = sr.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
						LOG.error("file is not sorted on coordinate :"+header.getSortOrder()+" "+path);
						return -1;
						}
					long count_bases = 0L;
					long sum_coverage = 0L;
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
								if(rec.getReadUnmappedFlag()) continue;
								if(rec.isSecondaryOrSupplementary()) continue;
								if(rec.getDuplicateReadFlag()) continue;
								if(rec.getMappingQuality() < this.mapping_quality ) continue;								
								}
							if(rec==null || !rec.getContig().equals(prevContig)) {
								if(coverage!=null) {//DUMP
									long count_bases_ctg = 0L;
									long sum_coverage_ctg = 0L;
									
									for(int i=0;i< coverage.length;i++) {
										if(mask.get(i)) continue;
										count_bases_ctg++;
										sum_coverage_ctg += coverage[i];
										}
									out.print(path);
									out.print("\t");
									out.print(sample);
									out.print("\t");
									out.print(prevContig);
									out.print("\t");
									out.print(count_bases_ctg);
									out.print("\t");
									out.print(sum_coverage_ctg);
									out.print("\t");
									if(count_bases_ctg>0) {
										out.print(sum_coverage_ctg/(double)count_bases_ctg);
										}
									else
										{
										out.print("N/A");
										}
									out.println();
									
									count_bases += count_bases_ctg;
									sum_coverage += sum_coverage_ctg;
									
									
									}
								coverage=null;
								mask=null;
								///
								System.gc();
								if(rec==null) break;
								
								final SAMSequenceRecord ssr = dict.getSequence(rec.getContig());
								coverage = new int[ssr.getSequenceLength()];
								mask = new BitSet(ssr.getSequenceLength());
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
							
							for(final AlignmentBlock block:rec.getAlignmentBlocks()) {
								int pos1=block.getReferenceStart();
								final int len = block.getLength();
								for(int i=0;i< len;i++) {
									if(pos1>0 && pos1 <= coverage.length) {
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
					out.print(count_bases);
					out.print("\t");
					out.print(sum_coverage);
					out.print("\t");
					if(count_bases>0) {
						out.print(sum_coverage/(double)count_bases);
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
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
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
