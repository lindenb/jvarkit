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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.pcr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
BEGIN_DOC

## Motivation

Mark PCR reads to their PCR amplicon https://www.biostars.org/p/149687/

Experimental, use with care.

* reads must have a read group id (see picard  addorreplacereadgroups http://broadinstitute.github.io/picard/ )
* if the read is not a paired-end read , then the quality of the read is set to zero
* if the MATE is not mapped , then the quality of the current read is set to zero
* if the properly-paired flag is not set, then the quality of the current read is set to zero
* if the read and the mate are not mapped on the same chromosome, then the quality of the current read is set to zero
* if the read and the mate are mapped on the same strand, then the quality of the current read is set to zero
* if the read(POS) < mate(POS) and read is mapped on the negative strand, then the quality of the current read is set to zero
* if the read(POS) > mate(POS) and read is mapped on the positive strand, then the quality of the current read is set to zero
* if no BED fragment is found overlapping the read, then the quality of the read is set to zero
* if a PCR fragment is found, the read group 'ID' is changed to 'ID_fragmentname'


## Example

create a bed file with bedtools, add a column 'PCR%%' as the name of the PCR fragment. Add a read group to ex1.bam , change the readgroup with pcrslicereads, select mapped reads having a MAPQ>0, sort and index, get the coverage with GATK/DepthOfCoverage using **--partitionType readgroup **


```makefile
.SHELL=/bin/bash


coverage : clipped.bam
	java -jar gatk/3.3.0/GenomeAnalysisTK.jar \
	   -T DepthOfCoverage   -omitBaseOutput \
	   -R samtools-0.1.19/examples/ex1.fa \
	   --partitionType readgroup -I $< -o coverage
	head  coverage.read_group_summary
	tail coverage.read_group_summary

clipped.bam : dist/pcrslicereads.jar test.bed
	java -jar picard/picard-tools-1.126/picard.jar  AddOrReplaceReadGroups \
		I=samtools-0.1.19/examples/ex1.bam \
		O=test.bam ID=myid LB=mylb PL=illumina PU=mypu SM=mysample
	java -jar $< -c   -B  test.bed -b test.bam |\
	samtools  view -q 1 -F 4 -bu  -  |\
	samtools  sort - clipped && samtools index clipped.bam

dist/pcrslicereads.jar:   src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrSliceReads.java
	make pcrslicereads

test.bed: samtools-0.1.19/examples/ex1.fa.fai
	bedtools makewindows -g $< -w 200 -s 50 |\
		awk '{printf("%s\tPCR%d\n",$$0,NR);}' > $@


```

```
$ head test.bed 
seq1    0       200     PCR1
seq1    50      250     PCR2
seq1    100     300     PCR3
seq1    150     350     PCR4
seq1    200     400     PCR5
seq1    250     450     PCR6
seq1    300     500     PCR7
seq1    350     550     PCR8
seq1    400     600     PCR9
seq1    450     650     PCR10
```

```
$ samtools view -h clipped.bam | more
@HD     VN:1.4  SO:coordinate
@SQ     SN:seq1 LN:1575 UR:file:/commun/data/packages/samtools-0.1.19/examples/ex1.fa   AS:ex1  M5:426e31835a6dfdcbf6c534671edf02f7
@SQ     SN:seq2 LN:1584 UR:file:/commun/data/packages/samtools-0.1.19/examples/ex1.fa   AS:ex1  M5:b6853ffe730ece50076db834dea18e3b
@RG     ID:myid PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR1    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR2    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR3    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR4    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR5    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR6    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR7    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR8    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR9    PL:illumina     PU:mypu LB:mylb SM:mysample
(...)
EAS54_71:8:105:854:975  83      seq2    1523    71      28M5S   =       1354    -202    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTG       <<<<;<:<<;<&<;<<<<<<<<<<<<<<<<<<<       H0:i:85 H1:i:85 MF:i:18 R
G:Z:myid_PCR60  NM:i:0  UQ:i:0  Aq:i:0
EAS139_11:7:50:1229:1313        83      seq2    1528    77      23M12S  =       1376    -187    TTTTTTCTTTTTTTTTTTTTTTTTTTTGCATGCCA     <<<<,<&<7<<<<<<<<<<<<<<<<<<<<<<<<<<     H0:i:3  H1:i:7  M
F:i:18  RG:Z:myid_PCR60 NM:i:1  UQ:i:11 Aq:i:0
EAS54_65:3:320:20:250   147     seq2    1532    77      19M16S  =       1367    -200    TTTTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAA     +'''/<<<<7:;+<;::<<<;;<<<<<<<<<<<<<     H0:i:1  H1:i:2  MF:i:18 R
G:Z:myid_PCR60  NM:i:2  UQ:i:24 Aq:i:6
EAS114_26:7:37:79:581   83      seq2    1533    68      18M17S  =       1349    -219    TTTTTTTTTTTTTTTTTTTTTTTCATGCCAGAAAA     3,,,===6===<===<;=====-============     H0:i:0  H1:i:1  MF:i:18 R
G:Z:myid_PCR60  NM:i:2  UQ:i:23 Aq:i:27
```

```
$ head coverage.read_group_summary
sample_id       total   mean    granular_third_quartile granular_median granular_first_quartile %_bases_above_15
mysample_rg_myid_PCR2   1106    0.35    1       1       1       0.2
mysample_rg_myid_PCR3   674     0.21    1       1       1       0.0
mysample_rg_myid_PCR1   51      0.02    1       1       1       0.0
mysample_rg_myid_PCR6   1279    0.40    1       1       1       0.1
mysample_rg_myid_PCR7   1554    0.49    1       1       1       1.3
mysample_rg_myid_PCR4   1147    0.36    1       1       1       0.8
mysample_rg_myid_PCR5   1738    0.55    1       1       1       1.5
mysample_rg_myid_PCR9   1260    0.40    1       1       1       0.8
mysample_rg_myid_PCR8   1930    0.61    1       1       1       2.1
```


END_DOC
 *
 */
@Program(
	name="pcrslicereads",
	description="Mark PCR reads to their PCR amplicon",
	keywords={"pcr","sam","bam","cigar"},
	biostars=149687
	)
public class PcrSliceReads extends Launcher
	{
	private static final Logger LOG =Logger.build(PcrSliceReads.class).make();
	
	@Parameter(names="-c",description="clip read to their PCR fragments. see https://github.com/lindenb/jvarkit/wiki/PcrClipReads")
	private boolean clip_reads=false;
	private IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	private enum AmbiguityStrategy {zero,random,closer};
	@Parameter(names="-a",description=" if a read is mapped on multiple PCR fragments, how to resolve ambiguity")
	private AmbiguityStrategy ambiguityStrategy= AmbiguityStrategy.closer; 
	@Parameter(names="--random",description=" random seed")
	private Random random=RandomConverter.now();//0L for reproductive calculations
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-B","--bed"},description="bed file containing non-overlapping PCR fragments. Column name is required.")
	private File bedFile=null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	private List<Interval> findIntervals(String chrom,int start,int end)
		{
		Interval i= new Interval(chrom,start,end);
		return  new ArrayList<Interval>(this.bedIntervals.getOverlapping(i));
		}
	
	
	@SuppressWarnings("resource")
	private int run(SamReader reader)
		{
		ReadClipper readClipper = new ReadClipper();
		SAMFileHeader header1= reader.getFileHeader();
		SAMFileHeader header2 = header1.clone();
		header2.addComment(getProgramName()+" "+getVersion()+": Processed with "+getProgramCommandLine());
		header2.setSortOrder(SortOrder.unsorted);
		
		for(SAMReadGroupRecord srg:header1.getReadGroups())
			{
			
			for(Interval i:this.bedIntervals.keySet())
				{
				//create new read group for this interval
				SAMReadGroupRecord rg=new SAMReadGroupRecord(srg.getId()+"_"+this.bedIntervals.get(i).getName(), srg);
				header2.addReadGroup(rg);
				}
			}
		
		SAMFileWriter sw=null;
		SAMRecordIterator iter = null;
		try
			{
			sw=writingBamArgs.openSAMFileWriter(outputFile, header2, false);
			SAMSequenceDictionaryProgress progress =new SAMSequenceDictionaryProgress(header1);
			iter =  reader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec= progress.watch(iter.next());
				if(rec.getReadUnmappedFlag())
					{
					sw.addAlignment(rec);
					continue;
					}
				
				if(!rec.getReadPairedFlag())
					{
					//@doc if the read is not a paired-end read ,  then the quality of the read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getMateUnmappedFlag())
					{
					//@doc if the MATE is not mapped ,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				if(!rec.getProperPairFlag())
					{
					//@doc if the properly-paired flag is not set,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getMateReferenceIndex()!=rec.getReferenceIndex())
					{
					//@doc if the read and the mate are not mapped on the same chromosome,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag())
					{
					//@doc if the read and the mate are mapped on the same strand,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				int chromStart;
				int chromEnd;
				if(rec.getAlignmentStart() < rec.getMateAlignmentStart())
					{
					if(rec.getReadNegativeStrandFlag())
						{
						//@doc if the read(POS) < mate(POS) and read is mapped on the negative strand, then the quality of the current read is set to zero
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					else
						{
						chromStart = rec.getAlignmentStart();
						chromEnd = rec.getMateAlignmentStart();
						}
					}
				else 
					{
					if(!rec.getReadNegativeStrandFlag())
						{
						//@doc if the read(POS) > mate(POS) and read is mapped on the positive strand, then the quality of the current read is set to zero
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					else
						{
						chromStart = rec.getMateAlignmentStart();
						chromEnd = rec.getAlignmentStart();//don't use getAlignmentEnd, to be consistent with mate data
						}
					}
				
				
				
				
				List<Interval> fragments = findIntervals(rec.getContig(),chromStart,chromEnd);
				if(fragments.isEmpty())
					{
					//@doc if no BED fragment is found overlapping the read, then the quality of the read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				Interval best=null;
				if(fragments.size()==1)
					{
					best = fragments.get(0);
					}
				else switch(this.ambiguityStrategy)
					{
					case random:
						{
						best = fragments.get(this.random.nextInt(fragments.size()));
						break;
						}
					case zero:
						{
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					case closer:
						{
						int best_distance=Integer.MAX_VALUE;
						for(Interval frag : fragments)
							{
							int distance= Math.abs(chromStart-frag.getStart()) + Math.abs(frag.getEnd() - chromEnd);
							if(best==null || distance < best_distance)
								{
								best = frag;
								best_distance = distance;
								}
							}
						break;
						}
					default: throw new IllegalStateException();
					}
				
				
				
				if(clip_reads)
					{
					rec = readClipper.clip(rec, best);
					if(rec.getMappingQuality()==0)
						{
						sw.addAlignment(rec);
						continue;
						}
					}
				SAMReadGroupRecord rg = rec.getReadGroup();
				if(rg == null )
					{
					throw new IOException("Read "+rec+" is missing a Read-Group ID . See http://broadinstitute.github.io/picard/ http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups");
					}
				rec.setAttribute("RG",rg.getId()+"_"+best.getName());
				sw.addAlignment(rec);
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sw);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
	if(bedFile==null)
			{
			LOG.error("undefined bed file");
			return -1;
			}
		BufferedReader r=null;
		SamReader samReader=null;
		try {
			samReader = super.openSamReader(oneFileOrNull(args));
			
			final BedLineCodec codec=new BedLineCodec();
			 r= IOUtils.openFileForBufferedReading(bedFile);
			String line;
			while((line=r.readLine())!=null)
				{
				final BedLine bed=codec.decode(line);
				if(bed==null) continue;
				final String chrom =  bed.getContig();
				int chromStart1 = bed.getStart();
				int chromEnd1 =  bed.getEnd();
				if(chromStart1<1 || chromStart1>chromEnd1)
					{
					LOG.error("Bad bed line "+line);
					return -1;
					}
				final String name = bed.get(3).trim();
				if(name==null || name.isEmpty())
					{
					LOG.error("Bad bed line (name missing) in  "+line);
					return -1;
					}
				final Interval i =new Interval(chrom, chromStart1, chromEnd1,false,name);
				this.bedIntervals.put(i, i);
				}
			return run(samReader);
			}
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(samReader);
			}
		}

	
	public static void main(String[] args) {
		new PcrSliceReads().instanceMainWithExit(args);
		}

}
