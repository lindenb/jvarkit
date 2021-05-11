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
package com.github.lindenb.jvarkit.tools.bamxdepth;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.IntPredicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
/**
BEGIN_DOC

## Input

input is a set of indexed BAM/CRAM files or a file with the '.list' suffix containing the paths to the BAM/CRAM paths.

## Example

```
$ java -jar dist/bamxtremdepth.jar --interval_list -M 15 -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null 
@HD	VN:1.6
@SQ	SN:RF01	LN:3302	M5:59dccb944425dd61f895a564ad7b56a7	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF02	LN:2687	M5:2c9c1ac1b7468ffaae96ad91c095c8b5	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF03	LN:2592	M5:d41f980f20d9cefbfd11ba2c1078f078	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF04	LN:2362	M5:935a2ad8d2f573d50b8c427f3b2f7c5d	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF05	LN:1579	M5:cdbaf0c352b44a79bef98daba1940d8a	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF06	LN:1356	M5:bdad6ae494bde13504d6988f2c7d94cd	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF07	LN:1074	M5:fa90eeb699fb6090b8f94611f19ac219	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF08	LN:1059	M5:4f8b3f714dc2327655941e6386c96b3b	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF09	LN:1062	M5:b221800f99aa2dc42258c0011ec228c8	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF10	LN:751	M5:f504f07ea3564b984207376aa02d8d00	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF11	LN:666	M5:7a7cf2c7813f2e8bd74be383014202ca	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@CO	bamxtremdepth. compilation:20210511100604 githash:47c5dc0 htsjdk:2.23.0 date:20210511101013. cmd:--interval_list -M 15 -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
RF02	1	6	+	DP<=0
RF02	2686	2687	+	DP<=0
RF03	1	1	+	DP<=0
RF04	1	9	+	DP<=0
RF04	2358	2362	+	DP<=0
RF05	1	1	+	DP<=0
RF05	543	543	+	DP>=15
RF05	1573	1579	+	DP<=0
RF06	1	2	+	DP<=0
RF07	1	5	+	DP<=0
RF07	1074	1074	+	DP<=0
RF09	1	1	+	DP<=0
RF09	1062	1062	+	DP<=0
RF10	1	1	+	DP<=0
RF10	362	373	+	DP<=0
RF10	751	751	+	DP<=0
RF11	64	72	+	DP>=15
RF11	74	74	+	DP>=15
RF11	296	304	+	DP<=0
RF11	590	590	+	DP>=15
RF11	595	596	+	DP>=15

```


END_DOC
 */

@Program(name="bamxtremdepth",
description="Compute low and high depth shared by a set of BAM files",
keywords={"bam","depth","coverage"},
creationDate="20210511",
modificationDate="20210511",
generate_doc=true
)
public class BamXtremDepth extends Launcher {
	private static final Logger LOG = Logger.build(BamXtremDepth.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxPath =null;
	@Parameter(names={"--validation-stringency"},description="SAM Reader Validation Stringency")
	private ValidationStringency validationStringency = ValidationStringency.LENIENT;
	@Parameter(names={"-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile =null;
	@Parameter(names={"-Q","--mapq"},description="min mapping quality")
	private int mapq=1;
	@Parameter(names={"--min-coverage","-m"},description="inclusive min coverage")
	private int min_coverage= 0;
	@Parameter(names={"--max-coverage","-M"},description="inclusive max coverage")
	private int max_coverage= 1_000;
	@Parameter(names={"--interval_list"},description="force interval list output (default is BED)")
	private boolean output_format_interval_list = false;
	@Parameter(names={"--contig"},description="limit to this chromosome.")
	private Set<String> limit_to_chromosomes = new HashSet<>();
	@Parameter(names={"-A"},description="alternative algorithm. Scan the whole chromosome instead of scanning the previously found interval. Could be faster if there are too many intervals.")
	private boolean alternative_algorithm = false; 	

	/**
	 * @param srcLoc the previous intervals
	 * @param coverage the coverage calculated on srcLoc
	 * @param depthTest test the coverage
	 * @return intersection of srcLoc and (coverage+depthTest)
	 */
	private void reduce(
			final BitSet bitSet,
			final CoverageFactory.SimpleCoverage coverage,
			final IntPredicate depthTest) {
		for(int x= coverage.getStart(); x<= coverage.getEnd();++x) {
			if(!bitSet.get(x-1)) continue;
			int depth = coverage.get(x-coverage.getStart());
			if(!depthTest.test(depth)) bitSet.set(x-1,false);
			}
		}
	
	
	private void runLoop(
			final List<Path> bamPaths,
			final SamReaderFactory srf,
			final SAMSequenceRecord ssr,
			final SAMSequenceDictionary refDict,
			BitSet lowList,
			BitSet  highList
			) throws IOException{
		final CoverageFactory coverageFactory = new CoverageFactory().setMappingQuality(mapq);
		for(final Path bamPath: bamPaths) {		
			try(SamReader samReader =srf.open(bamPath)) {
				// check indexed
				if(!samReader.hasIndex()) {
					throw new JvarkitException.BamHasIndex(bamPath.toString());
					}
				final SAMFileHeader header = samReader.getFileHeader();
				// check dictionary
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				SequenceUtil.assertSequenceDictionariesEqual(dict, refDict);
				// get the coverage in the low and high region.
				final CoverageFactory.SimpleCoverage coverage;
				if (this.alternative_algorithm) {
					LOG.info(ssr.getContig()+" "+bamPath);
					coverage = coverageFactory.getSimpleCoverage(samReader, ssr, null);
					}
				else	{
					final List<Locatable> regions = bitSetToLocatables(ssr, lowList);
                        		regions.addAll(bitSetToLocatables(ssr,highList));
                        		if(regions.isEmpty()) break;
					LOG.info(ssr.getContig()+" "+bamPath+ " n-intervals:"+regions.size() + " size:"+regions.stream().mapToInt(R->R.getLengthOnReference()).sum());
				 	coverage = coverageFactory.getSimpleCoverage(samReader, regions, null);
					}
				// compute the new low and hight coverage
				reduce(lowList, coverage, I->I<=this.min_coverage);
				reduce(highList, coverage, I->I>=this.max_coverage);
				}
			System.gc();
			}
		}
	
	private List<Locatable> bitSetToLocatables(
		final SAMSequenceRecord ssr,
		final BitSet bitSet) {
		final List<Locatable> L = new ArrayList<>();
		int i=0;
		while(i< bitSet.size()) {
			i = bitSet.nextSetBit(i);
			if(i==-1) break;
			int j=i+1;
			while(j< bitSet.size() && bitSet.get(j)) {
				j++;
				}
			L.add(new SimpleInterval(ssr.getContig(), i+1, j));//bitSet is zero based
			i=j;
			}
		return L;
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		if(this.min_coverage>= this.max_coverage) {
			LOG.equals("min coverage >= max_coverage");
			return -1;
			}
		
		try {
			final SamReaderFactory srf  = SamReaderFactory.
				make().
				validationStringency(this.validationStringency).
				referenceSequence(this.faidxPath);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(faidxPath);
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			if(bamPaths.isEmpty()) {
				LOG.error("No bam was defined");
				return -1;
				}
			try (PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				final boolean out_is_intervallist = (
						this.output_format_interval_list || 
						(this.outputFile!=null && (this.outputFile.toString().endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST) || this.outputFile.toString().endsWith(FileExtensions.INTERVAL_LIST)))
						);
				if(out_is_intervallist) {
					final SAMFileHeader outHeader = new SAMFileHeader(dict);
					JVarkitVersion.getInstance().addMetaData(this, outHeader);
					final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
		            codec.encode(out, outHeader);
					}
			
				for(final SAMSequenceRecord ssr : dict.getSequences()) {
					if(!this.limit_to_chromosomes.isEmpty()  && !this.limit_to_chromosomes.contains(ssr.getSequenceName())) continue;
					LOG.info("Scanning "+ ssr.getSequenceName());
					final BitSet highList = new BitSet(ssr.getLengthOnReference());
					highList.set(0, ssr.getSequenceLength());
					final BitSet lowList = new BitSet(ssr.getLengthOnReference());
					lowList.set(0, ssr.getSequenceLength());
					
					runLoop(bamPaths, srf, ssr, dict, lowList, highList);
					
					final List<Interval> rgn = new ArrayList<>(highList.size()+lowList.size());
					bitSetToLocatables(ssr,lowList).stream().map(R->new Interval(R.getContig(),R.getStart(),R.getEnd(),false,"LE."+this.min_coverage)).forEach(R->rgn.add(R));
					bitSetToLocatables(ssr,highList).stream().map(R->new Interval(R.getContig(),R.getStart(),R.getEnd(),false,"GE."+this.max_coverage)).forEach(R->rgn.add(R));
					rgn.sort(new ContigDictComparator(dict).createLocatableComparator());
					for(final Interval R: rgn) {
						out.print(R.getContig());
						out.print("\t");
						out.print(R.getStart()-(out_is_intervallist?0:1));
						out.print("\t");
						out.print(R.getEnd());
						out.print("\t");
						if(out_is_intervallist) out.print("+\t");
						out.print(R.getName());
						out.println();
						}
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
		
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamXtremDepth().instanceMainWithExit(args);

		}

	}
