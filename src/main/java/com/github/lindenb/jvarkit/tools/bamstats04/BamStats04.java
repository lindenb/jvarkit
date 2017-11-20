/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.IntUnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
BEGIN_DOC

## History

* 2017-11-20: added new column 'partition'
* 2017-11-20: can read more than one BAM File.

## Example

```
$ java -jar dist/bamstats04.jar -B src/test/resources/toy.bed.gz src/test/resources/toy.bam 2> /dev/null | column -t 

#chrom  start  end  length  sample  mincov  maxcov  meancov  mediancov  nocoveragebp  percentcovered
ref     10     13   3       S1      3       3       3.0      3.0        0             100
ref2    1      2    1       S1      2       2       2.0      2.0        0             100
ref2    13     14   1       S1      6       6       6.0      6.0        0             100
ref2    16     17   1       S1      6       6       6.0      6.0        0             100
```

END_DOC
 */
@Program(name="bamstats04",
	description="Coverage statistics for a BED file.",
	keywords={"sam","bam","coverage","depth","statistics","bed"}
	)
public class BamStats04 extends Launcher
	{
	private static final Logger LOG = Logger.build(BamStats04.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-cov","--cov"},description="min coverage to say the position is not covered")
	private int MIN_COVERAGE = 0 ;
	@Parameter(names={"-f","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();
	@Parameter(names={"-B","--bed"},description="Bed File. Required",required=true)
	private File bedFile = null;
	@Parameter(names={"-R","--ref"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION+" If set, a column with the GC% will be added")
	private File faidxFile = null;
	@Parameter(names={"-partition","--partition"},description="[20171120]"+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	
	
	/** map depth to 0 if depth <= MIN_COVERAGE */
	private final IntUnaryOperator depthAdjuster = (D)->(D<=this.MIN_COVERAGE?0:D);
	
	private static class IntervalStat
		{	
		private final BedLine bedLine;
		private final int counts[];
		IntervalStat(final BedLine bedLine) {
			this.bedLine = bedLine;
			this.counts=new int[bedLine.getEnd()-bedLine.getStart()+1];
			
			Arrays.fill(this.counts, 0);
			}
		void visit(final SAMRecord rec) {
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
			int refpos1=rec.getAlignmentStart();
    		for(final CigarElement ce:cigar)
    			{
    			final CigarOperator op=ce.getOperator();
    			if(!op.consumesReferenceBases()) continue;
    			if(op.consumesReadBases())
    				{
    				for(int i=0;i< ce.getLength();++i)
		    			{
						if(refpos1+i>= this.bedLine.getStart() && refpos1+i<= this.bedLine.getEnd())
							{
							this.counts[refpos1+i-this.bedLine.getStart()]++;
							}
	    				}
    				}
    			refpos1+=ce.getLength();
    			if(refpos1>this.bedLine.getEnd()) break;
    			}
			}
		
		}
	
	@Override
		public int doWork(final List<String> args) {
			if(this.bedFile==null || !this.bedFile.exists()) {
				LOG.error("undefined option -B (bed file)");
				return -1;
			}
			if(args.isEmpty())  {
				LOG.error("Bam files missing");
				return -1;
				}
			final String NO_PARTITION="N/A";
			BufferedReader bedIn=null;
			final List<SamReader> samReaders = new ArrayList<>(args.size());
			PrintWriter pw = null;
			IndexedFastaSequenceFile indexedFastaSequenceFile=null;
			GenomicSequence genomicSequence=null;
			try
				{
				final BedLineCodec codec= new BedLineCodec();
				
				bedIn=IOUtils.openFileForBufferedReading(this.bedFile);
				SAMSequenceDictionary dict = null;
				
				for(final String filename: IOUtils.unrollFiles(args)) {
					LOG.info(filename);
					final SamReader samReader = super.openSamReader(filename);
					if(!samReader.hasIndex()) {
						LOG.error(filename+" is not indexed");
						samReader.close();
						return -1;
						}
					final SAMSequenceDictionary d = samReader.getFileHeader().getSequenceDictionary();
						if(d==null) {
						samReader.close();
						LOG.error("SAM sequence dictionary missing in SAM Header");
						return -1;
						}
					
					samReaders.add(samReader);
						
					if(dict==null) {
						dict=d;
						}
					else if(SequenceUtil.areSequenceDictionariesEqual(d, dict)) {
						throw new JvarkitException.DictionariesAreNotTheSame(d, dict);
						}
					}
				
				if(samReaders.isEmpty()) {
					LOG.error("No Bam defined");
					return -1;
				}
				
				
				
				
				if(this.faidxFile!=null) {
					indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidxFile);
				}
				pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
				pw.println("#chrom\tstart\tend\tlength\t"+
					this.partition.name()+"\t"+
					(indexedFastaSequenceFile==null?"":"gc_percent\t")+
					"mincov\tmaxcov\tmeancov\tmediancov\tnocoveragebp\tpercentcovered");
	
			
				String line=null;
				while((line=bedIn.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					final BedLine bedLine = codec.decode(line);
					if(bedLine==null) continue;
					if(dict.getSequence(bedLine.getContig())==null)
						{
						LOG.error("Unknown contig in "+line);
						return -1;
						}
					
					if(bedLine.getStart()>bedLine.getEnd())
						{
						LOG.info("ignoring "+bedLine);
						continue;
						}
					if(indexedFastaSequenceFile!=null && (genomicSequence==null || !genomicSequence.getChrom().equals(bedLine.getContig()))) {
						genomicSequence = new GenomicSequence(indexedFastaSequenceFile, bedLine.getContig());
						}
					
					final Map<String, IntervalStat> sample2stats= new HashMap<>();
					
					for(final SamReader samReader:samReaders) 
						{
						/**
						 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
			    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
						 */
						final SAMRecordIterator r=samReader.queryOverlapping(
								bedLine.getContig(),
								bedLine.getStart(),
								bedLine.getEnd()
								);
						while(r.hasNext())
							{
							final SAMRecord rec=r.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(this.filter.filterOut(rec)) continue;
							if(!rec.getReferenceName().equals(bedLine.getContig())) continue;
							
							
							final String partition;
							final SAMReadGroupRecord group = rec.getReadGroup();
							if(group==null)
								{
								partition=NO_PARTITION;
								}
							else
								{
								final String name = this.partition.apply(group);
								partition = (StringUtil.isBlank(name)?NO_PARTITION:name);
								}
							
							IntervalStat stat= sample2stats.get(partition);
							if(stat==null) 
								{
								stat = new IntervalStat(bedLine);
								sample2stats.put(partition,stat);
								}
							stat.visit(rec);
							}
						
						r.close();
						} // end of loop over sam Readers
					
					final Integer gcPercent = (genomicSequence==null?
						null:
						genomicSequence.getGCPercent(bedLine.getStart()-1,bedLine.getEnd())).getGCPercentAsInteger()
						;
					
					
					for(final String partitionName : sample2stats.keySet()) {
						final IntervalStat stat = sample2stats.get(partitionName);
						Arrays.sort(stat.counts);
						
						final int count_no_coverage=(int)Arrays.stream(stat.counts).
								filter(D->this.depthAdjuster.applyAsInt(D)<=0).
								count()
								;
						
						final double mean= Percentile.average().evaluate(Arrays.stream(stat.counts).
								map(this.depthAdjuster)
								);
						
		                final double median_depth = Percentile.median().evaluate(Arrays.stream(stat.counts).
								map(this.depthAdjuster)
								);
		                
						
						pw.println(
								bedLine.getContig()+"\t"+
								(bedLine.getStart()-1)+"\t"+
								(bedLine.getEnd())+"\t"+
								stat.counts.length+"\t"+
								partitionName+"\t"+
								(gcPercent==null? "":String.valueOf(gcPercent)+"\t")+
								stat.counts[0]+"\t"+
								stat.counts[stat.counts.length-1]+"\t"+
								mean+"\t"+
								median_depth+"\t"+
								count_no_coverage+"\t"+
								(int)(((stat.counts.length-count_no_coverage)/(double)stat.counts.length)*100.0)
								);
						
						}
					}
				pw.flush();
				pw.close();pw=null;
				LOG.info("done");
				return RETURN_OK;
				}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(pw);
			CloserUtil.close(bedIn);
			CloserUtil.close(samReaders);
			}
		}
	
	public static void main(final String[] args) throws Exception
		{
		new BamStats04().instanceMainWithExit(args);
		}
}

