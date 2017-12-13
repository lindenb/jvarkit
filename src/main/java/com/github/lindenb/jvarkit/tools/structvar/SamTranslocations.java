/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIterator;
import com.github.lindenb.jvarkit.util.iterator.ForwardPeekIteratorImpl;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
/**

BEGIN_DOC

## History

* 2017-12-13 :  refactoring for balanced translocation.

END_DOC
*/
@SuppressWarnings("unused")
@Program(name="samtranslocations",
	description="Explore balanced translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam","xslt","xml"}
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.build(SamTranslocations.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--region","--interval"},description="Limit analysis to this interval. "+IntervalParser.OPT_DESC)
	private String region_str=null;
	@Parameter(names={"--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamRecordJEXLFilter.buildAcceptAll();
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;

	@Parameter(names={"-md","--max-distance"},description="Max distance between forward-reverse")
	private int max_distance = 1000;

	private static class PartitionState
		{
		final String name;
		SAMRecord last_rec = null;
		PartitionState(final String name) {
			this.name=name;
			}
		}	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_distance<=0) {
			LOG.error("max_distance <=0 ("+max_distance+")");
			return -1;
		}
		PrintWriter pw = null;
		ConcatSam.ConcatSamIterator samIter = null;
		ForwardPeekIterator<SAMRecord> forwardPeekIterator = null;
		try {
			
			samIter = new ConcatSam.Factory().setInterval(this.region_str).open(args);
			final SAMSequenceDictionary dict = samIter.getFileHeader().getSequenceDictionary();
			if(dict.size()<2) {
				LOG.error("Not enough contigs in sequence dictionary. Expected at least 2.");
				return -1;
			}
			forwardPeekIterator = new ForwardPeekIteratorImpl<>(
					new FilterIterator<>(samIter,rec->{
						if(rec.getReadUnmappedFlag()) return false;
						if(!rec.getReadPairedFlag()) return false;
						if(rec.getMateUnmappedFlag()) return false;
						final int tid1 =  rec.getReferenceIndex();
						final int tid2 =  rec.getMateReferenceIndex();
						if(tid1==tid2) return false;
						if(this.samRecordFilter.filterOut(rec)) return false;
						return true;
						})
					);
			
			pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			pw.print("#chrom1");
			pw.print('\t');
			pw.print("chrom1-start");
			pw.print('\t');
			pw.print("chrom1-end");
			pw.print('\t');
			pw.print("middle1\tstrand_plus_before_mid1_count\tstrand_plus_before_mid1_percent\tstrand_minus_after_mid1_count\tstrand_minus_after_mid1_percent");
			pw.print('\t');
			pw.print("chrom2");
			pw.print('\t');
			pw.print("chrom2-start");
			pw.print('\t');
			pw.print("chrom2-end");
			pw.print('\t');
			pw.print("middle2\tstrand_plus_before_mid2_count\tstrand_plus_before_mid2_percent\tstrand_minus_after_mid2_count\tstrand_minus_after_mid2_percent");
			pw.print('\t');
			pw.print("count-reads");
			pw.print('\t');
			pw.print("count-clipped");
			pw.print('\t');
			pw.print("partition");
			pw.println();
			
			final Map<String,PartitionState> partition2state = new HashMap<>();
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(dict).logger(LOG);
			while(forwardPeekIterator.hasNext()) {
				final SAMRecord rec = progress.watch(forwardPeekIterator.next());
				if(rec.getReadNegativeStrandFlag()) continue;// searching for -->
				final String partition = this.samRecordPartition.getPartion(rec, "N/A");
				PartitionState partitionState = partition2state.get(partition);
				if(partitionState==null) {
					partitionState = new PartitionState(partition);
					partition2state.put(partition, partitionState);
					}
				if(partitionState.last_rec!=null)
					{
					if(partitionState.last_rec==rec) {
						//reset last
						partitionState.last_rec=null;
						}
					continue;
					}
				final  List<SAMRecord> recordList=new ArrayList<>();
				recordList.add(rec);
				
				int x=0;
				for(;;)
					{
					final SAMRecord rec2 = forwardPeekIterator.peek(x);
					if(rec2==null) {
						break;
					}
					if(!rec.getReferenceIndex().equals(rec2.getReferenceIndex())) {
						break;
						}
					if(rec2.getStart()-rec.getEnd()>this.max_distance) {
						break;
						}
					if(!partition.equals(this.samRecordPartition.getPartion(rec2, "N/A"))) {
						++x;						
						continue;
						}
					
					if(!rec2.getMateReferenceName().equals(rec.getMateReferenceName())) {
						++x;
						continue;
						}
					if(Math.abs(rec2.getMateAlignmentStart()-rec.getMateAlignmentStart())>this.max_distance) {
						++x;
						continue;
						}
					
					recordList.add(rec2);
					partitionState.last_rec = rec2;
					++x;
					}
				
				if(recordList.size()<2) continue;
								
				// chrom start stuff
				final List<SAMRecord> listF1=  recordList.stream().
							filter(R->!R.getReadNegativeStrandFlag()).
							collect(Collectors.toList());
				final List<SAMRecord> listR1=  recordList.stream().
						filter(R->R.getReadNegativeStrandFlag()).
						collect(Collectors.toList());
				if(listR1.isEmpty()) continue;
				
				
				final int start1 = (int)Percentile.median().evaluate(listF1.stream().mapToInt(R->R.getEnd()));
				final int end1 =   (int)Percentile.median().evaluate(listR1.stream().mapToInt(R->R.getStart()));
				final int mid1= (start1+end1)/2;
			
				// chrom end stuff
				final List<SAMRecord> listF2=  recordList.stream().
							filter(R->!R.getMateNegativeStrandFlag()).
							collect(Collectors.toList());
				final List<SAMRecord> listR2=  recordList.stream().
						filter(R->R.getMateNegativeStrandFlag()).
						collect(Collectors.toList());
				if(listF2.isEmpty()) continue;
				if(listR2.isEmpty()) continue;
				
				final int mid2= (int)Percentile.median().evaluate(listF2.stream().mapToInt(R->R.getMateAlignmentStart()));

				
				
				pw.print(rec.getReferenceName());
				pw.print('\t');
				pw.print(listF1.stream().mapToInt(R->R.getEnd()).max().getAsInt());
				pw.print('\t');
				pw.print(listR1.stream().mapToInt(R->R.getStart()).min().getAsInt());
				pw.print('\t');
				pw.print(mid1);
				pw.print('\t');
				pw.print(listF1.stream().filter(R->R.getEnd()<=mid1).count());
				pw.print('\t');
				pw.print((int)(100.0*(listF1.stream().filter(R->R.getEnd()<=mid1).count())/listF1.size()));
				pw.print('\t');
				pw.print(listR1.stream().filter(R->R.getStart()>=mid1).count());
				pw.print('\t');
				pw.print((int)(100.0*(listR1.stream().filter(R->R.getStart()>=mid1).count())/listR1.size()));
				pw.print('\t');
				
				pw.print(rec.getMateReferenceName());
				pw.print('\t');
				pw.print(recordList.stream().mapToInt(R->R.getMateAlignmentStart()).min().getAsInt());
				pw.print('\t');
				pw.print(recordList.stream().mapToInt(R->R.getMateAlignmentStart()).max().getAsInt());
				pw.print('\t');
				pw.print(mid2);
				pw.print('\t');
				pw.print(listF2.stream().filter(R->R.getMateAlignmentStart()<=mid2).count());
				pw.print('\t');
				pw.print((int)(100.0*(listF2.stream().filter(R->R.getMateAlignmentStart()<=mid2).count())/listF2.size()));
				pw.print('\t');
				pw.print(listR2.stream().filter(R->R.getMateAlignmentStart()>=mid2).count());
				pw.print('\t');
				pw.print((int)(100.0*(listR2.stream().filter(R->R.getMateAlignmentStart()>=mid2).count())/listR2.size()));
				pw.print('\t');
				pw.print(recordList.size());
				pw.print('\t');
				pw.print(recordList.stream().filter(R->R.getCigar()!=null && R.getCigar().isClipped()).count());
				pw.print('\t');
				pw.print(partition);
				pw.println();
				}
			progress.finish();
			forwardPeekIterator.close();forwardPeekIterator=null;
			samIter.close();samIter=null;
			pw.flush();
			pw.close();pw=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(forwardPeekIterator);
			CloserUtil.close(samIter);
			CloserUtil.close(pw);
			}
		}
	
	public static void main(String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
