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
package com.github.lindenb.jvarkit.tools.bam2wig;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

/**
BEGIN_DOC

## Input

input is stdin; or  one or more BAM file sorted on coordinate; or a file ending with '.list' and containing the PATH to some bam files.

About Wiggle: [https://genome.ucsc.edu/goldenpath/help/wiggle.html](https://genome.ucsc.edu/goldenpath/help/wiggle.html)

About BedGraph: [https://genome.ucsc.edu/goldenpath/help/bedgraph.html](https://genome.ucsc.edu/goldenpath/help/bedgraph.html)

## Memory

warning: the program is memory consuming, it allocates on array of integer of the size of your longest contig.

## History:

20171115: removed cast_to_integer replaced by 'format', added percentile. Removed options --zerolength and --mindepth.

## Aggregators:

* COVERAGE :  coverage, all sample merged
* CLIPPING : consider only clipped base
* INSERTION: consider only Cigar events I. Only one base in the reference is flagged
* DELETION: consider Cigar events 'N' and 'D':
* READ_GROUPS : Number of 'samples' having a depth greater than 'min-depth'
* CASE_CTRL: ratio median(coverage-cases)/median(coverage-controls)

## Example
the input file

```bash
java -jar dist/bam2wig.jar -w 1 -s 3   examples/toy.bam
```

```
track type=wiggle_0 name="__REPLACE_WIG_NAME__" description="__REPLACE_WIG_DESC__"
fixedStep chrom=ref start=7 step=3 span=1
1
3
3
3
1
1
0
0
1
0
2
2
1
fixedStep chrom=ref2 start=1 step=3 span=1
1
2
3
4
5
6
6
5
4
3
3
```

END_DOC
 */
@Program(name="bam2wig",
description="Bam to fixedStep Wiggle converter , or BED GRAPH. Parses the cigar String to get the depth. Memory intensive: must alloc sizeof(int)*size(chrom)",
keywords={"bam","wig","wiggle","bed"},
modificationDate="20190417"
)
public class Bam2Wig extends Launcher
	{
	private static final String UCSC_HEADER="track type=track_type name=\"__REPLACE_WIG_NAME__\" description=\"__REPLACE_WIG_DESC__\"";
	private static final Logger LOG = Logger.build(Bam2Wig.class).make();
	private enum WHAT {COVERAGE,CLIPPING,INSERTION,DELETION,READ_GROUPS,CASE_CTRL};

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-t","--header"},description="print a UCSC custom track header: something lile "+UCSC_HEADER+". Use `sed` to replace the tokens. e.g: `sed '/^track/s/__REPLACE_WIG_NAME__/My data/'` ")
	private boolean custom_track = false;
	@Parameter(names={"-s","--windowShift"},description="window shift")
	private int win_shift = 25 ;
	@Parameter(names={"-w","--windowSize"},description="window size")
	private int window_span = 100 ;
	@Parameter(names={"-f","--format"},description="`Printf` Format for the values. see https://docs.oracle.com/javase/tutorial/java/data/numberformat.html . Use \"%.01f\" to print an integer. \"%e\" for scientific notation.")
	private String printfFormat = "%.3f";
	@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"--percentile"},description="How to group data in the sliding window ?")
	private Percentile.Type percentilType = Percentile.Type.AVERAGE;
	@Parameter(names={"-bg","--bedgraph"},description="Produce a BED GRAPH instead of a WIGGLE file.")
	private boolean bedGraph= false;
	@Parameter(names={"--display"},description="What kind of data should we display ?")
	private WHAT whatDisplay= WHAT.COVERAGE;
	@Parameter(names={"--partition"},description="When using display READ_GROUPS, how should we partition the ReadGroup ? "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition= SAMRecordPartition.sample;
	@Parameter(names={"--mindepth","--mindp"},description="When using display READ_GROUPS, What is the minimal read depth that should be considered ?")
	private int min_depth=0;
	@Parameter(names={"--region","--interval"},description="Limit analysis to this interval. "+IntervalParser.OPT_DESC)
	private String region_str=null;
	@Parameter(names={"--pedigree","-ped"},description="Pedigree file for CASE_CTRL. " + Pedigree.OPT_DESCRIPTION )
	private File pedigreeFile=null;


	public Bam2Wig()
		{
		
		}
	
	private static abstract class Aggregator
		{
		abstract void visit(final int array[],final SAMRecord rec);
		void finish(final int array[]) {}
		protected void incr(final int array[],int pos1,int length)
			{
			for(int i=0;i< length ;++i)
				{
				final int array_index= pos1+i-1;
				if(array_index>0 && array_index<array.length)
					{
					array[array_index]++;
					}
				}
			}
		}
	
	private static class CoverageAggregator extends Aggregator
		{
		private final Predicate<SAMRecord> internalFilter;
		CoverageAggregator() {
			this(R->true);
			}
		protected CoverageAggregator(final Predicate<SAMRecord> predicate) {
			this.internalFilter = predicate;
			}
		@Override
		void visit(final int array[],final SAMRecord rec)
			{
			if(!this.internalFilter.test(rec)) return ;
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
    		int refpos1=rec.getAlignmentStart();
    		for(final CigarElement ce:cigar)
    			{
    			final CigarOperator op = ce.getOperator();
    			if(op.consumesReferenceBases())
    				{
    				if(op.consumesReadBases())
    					{
    					incr(array,refpos1,ce.getLength());
    					}
    				refpos1+=ce.getLength();
    				}    				
    			}
			}
		}
	
	private static class DeletionAggregator extends Aggregator
		{
		@Override
		void visit(final int array[],final SAMRecord rec)
			{
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
			int refpos1=rec.getAlignmentStart();
			for(final CigarElement ce:cigar)
				{
				final CigarOperator op = ce.getOperator();
				switch(op)
					{
					case D:
					case N:incr(array,refpos1,ce.getLength());
					default: break;
					}
				if(op.consumesReferenceBases())
					{
					refpos1+=ce.getLength();
					}
				}
			}
		}
	private static class InsertionAggregator extends Aggregator
		{
		@Override
		void visit(final int array[],final SAMRecord rec)
			{
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
			int refpos1=rec.getAlignmentStart();
			for(final CigarElement ce:cigar)
				{
				final CigarOperator op = ce.getOperator();
				switch(op)
					{
					case I: incr(array,refpos1,1);
					default: break;
					}
				if(op.consumesReferenceBases())
					{
					refpos1+=ce.getLength();
					}
				}
			}
		}
	
	
	private static class ClipAggregator extends Aggregator
		{
		@Override
		void visit(final int array[],final SAMRecord rec)
			{
			final Cigar cigar=rec.getCigar();
			if(cigar==null) return;
			int refpos1=rec.getUnclippedStart();//<- warning 
			for(final CigarElement ce:cigar)
				{
				final CigarOperator op = ce.getOperator();
				if(op.isClipping())
					{
					incr(array,refpos1,ce.getLength());
					refpos1+=ce.getLength();
					}
				else if(op.consumesReferenceBases())
					{
					refpos1+=ce.getLength();
					}
				}
			}
		}
	private static abstract class BufferedAggregator extends Aggregator
		{
		private final List<SAMRecord> buffer= new ArrayList<>();
		private int last_start=1;
		protected abstract void dump(final List<SAMRecord> records,int start1,int end1,final int array[]);
		
		protected abstract String partition(SAMRecord rec);
		
		protected Map<Integer,Counter<String>> fillPositions(final List<SAMRecord> records,final int begin1,final int end1)
			{
			final Map<Integer,Counter<String>> pos2sample2depth = new HashMap<>();
			for(final SAMRecord rec:records) {
				final Cigar cigar = rec.getCigar();
				if(cigar==null) continue;
				final String sample = partition(rec);
				if(StringUtil.isBlank(sample)) continue;
				
				int pos1= rec.getAlignmentStart();
				for(final CigarElement ce:cigar) {
					if(pos1>=end1) break;
					final CigarOperator op= ce.getOperator();
					if(op.consumesReferenceBases())
						{
						final int L=ce.getLength();
						if(op.consumesReadBases())
							{
							for(int i=0;i< L && pos1+i<end1;++i)
								{
								final int nx = pos1 + i;
								if(nx>=begin1 && nx< end1)
									{
									Counter<String> sample2coverage = pos2sample2depth.get(nx);
									if(sample2coverage==null) {
										sample2coverage = new Counter<String>();
										pos2sample2depth.put(nx, sample2coverage);
										}
									sample2coverage.incr(sample);
									}
								}
							}
						pos1+=L;
						}
					}
				
				}
			return pos2sample2depth;
			}
		
		@Override
		void visit(final int array[],final SAMRecord rec)
			{
			if(rec.getAlignmentStart() < this.last_start)
				{
				throw new IllegalStateException("got read "+rec+" after last_start="+last_start);
				}
			this.buffer.removeIf(R->R.getEnd()< this.last_start);
			final List<SAMRecord> toDump= buffer.stream().
					filter(R->R.getEnd()< R.getAlignmentEnd()).
					collect(Collectors.toList());
			
			if(!toDump.isEmpty()) {
				dump(toDump,this.last_start,rec.getAlignmentStart(),array);
				this.last_start = rec.getAlignmentStart();
				}
			this.buffer.add(rec);
			}
		@Override
		void finish(final int array[]) {
			dump(this.buffer,this.last_start,array.length+1,array);
			buffer.clear();
			this.last_start=1;
			}
		}
	
	private static class NumberOfSamplesCoveredX extends BufferedAggregator
		{
		private final int minDepth;
		private SAMRecordPartition samRecordPartition;
		NumberOfSamplesCoveredX(int minDepth, final SAMRecordPartition samRecordPartition) {
			this.minDepth = minDepth;
			this.samRecordPartition = samRecordPartition;
			}
		
		@Override
		protected String partition(final SAMRecord rec) {
			return this.samRecordPartition.apply(rec.getReadGroup());
			}
		
		@Override
		protected void dump(final List<SAMRecord> records,final int begin1,final int end1,final int array[]) {
			final Map<Integer,Counter<String>> pos2sample2depth= fillPositions(records,begin1,end1);
			for(int start1 = begin1;start1<end1;++start1)
				{
				final Counter<String> sample2coverage = pos2sample2depth.get(start1);
				if(sample2coverage==null) continue;
				
				final int num_samples= (int)sample2coverage.stream().
						mapToInt(E->E.getValue().intValue()).
						filter(D->D>=this.minDepth).
						count()
						;				
				if(start1>0 && start1<=array.length)
					{
					array[start1-1]=num_samples;
					}
				
				}
			}
		}
	
	private static class CaseControlAggregator extends BufferedAggregator
		{
		private final Map<String,Pedigree.Person> case2person;
		private final Map<String,Pedigree.Person> ctrl2person;
		CaseControlAggregator(final File pedigreeFile) {
			final Pedigree pedigree ;
			IOUtil.assertFileIsReadable(pedigreeFile);
			try {
				pedigree = Pedigree.newParser().parse(pedigreeFile);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			this.case2person=pedigree.getPersons().stream().
						filter(P->P.isAffected()).
						collect(Collectors.toMap(P->P.getId(), P->P))
						;
			this.ctrl2person=pedigree.getPersons().stream().
					filter(P->P.isUnaffected()).
					collect(Collectors.toMap(P->P.getId(), P->P))
					;
			}
		@Override
		void visit(int[] array, SAMRecord rec) {
			if(this.case2person.isEmpty()) return;
			if(this.ctrl2person.isEmpty()) return;
			if(StringUtil.isBlank(partition(rec))) return;
			super.visit(array, rec);
			}
		
		@Override
		protected String partition(final SAMRecord rec) {

			final SAMReadGroupRecord rg = rec.getReadGroup();
			if(rg==null) return null;
			final String sample= rg.getSample();			
			if(!(this.ctrl2person.containsKey(sample)|| this.case2person.containsKey(sample))) return null;
			return null;
			}
		
		@Override
		protected void dump(final List<SAMRecord> records,final int begin1,final int end1,final int array[]) {
			if(this.case2person.isEmpty()) return;
			if(this.ctrl2person.isEmpty()) return;
			final Map<Integer,Counter<String>> pos2sample2depth= fillPositions(records,begin1,end1);
				
					
			for(int start1 = begin1;start1<end1;++start1)
				{
				
				final Counter<String> sample2coverage = pos2sample2depth.get(start1);
				if(sample2coverage==null) continue;
				
				final double median_cases = Percentile.median().evaluate(this.case2person.keySet().
						stream().
						mapToDouble(ID->sample2coverage.count(ID))
						);
				final double median_ctrl = Percentile.median().evaluate(this.ctrl2person.keySet().
						stream().
						mapToDouble(ID->sample2coverage.count(ID))
						);

				final double ratio = median_cases / median_ctrl;
				
				if(start1>0 && start1<=array.length)
					{
					array[start1-1]= (int)(ratio * 1000.0);
					}
				
				}
			}
		}

	
	
	private void run(
			final PrintWriter pw,
			final CloseableIterator<SAMRecord> iter,
			final SAMSequenceDictionary dict,
			final Interval interval // may be null
			)
		{
		final Aggregator aggregator;
		switch(this.whatDisplay)
			{
			case COVERAGE: aggregator = new CoverageAggregator();break;
			case CLIPPING : aggregator = new ClipAggregator(); break;
			case INSERTION : aggregator = new InsertionAggregator();break;
			case DELETION : aggregator = new DeletionAggregator();break;
			case READ_GROUPS: aggregator = new NumberOfSamplesCoveredX(this.min_depth, this.partition);break;
			case CASE_CTRL : 
				if(this.pedigreeFile==null) {
					throw new JvarkitException.UserError("undefined pedigree");
				}
				aggregator = new CaseControlAggregator(this.pedigreeFile); break;
			default: throw new IllegalStateException(this.whatDisplay.name());
			}
		
		final Percentile percentile = Percentile.of(this.percentilType);
		SAMSequenceRecord ssr = null;
		int array[]=null;
		final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(dict);
		if(this.custom_track)
			{
			pw.println(
				UCSC_HEADER.replace("track_type", 
					this.bedGraph?"bedGraph":"wiggle_0")
					);
			}
		
		for(;;)
			{
			final SAMRecord rec;
			if(iter.hasNext())
				{
				rec=progess.watch(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(this.samRecordFilter.filterOut(rec)) continue;
				if(interval!=null && !interval.overlaps(rec)) continue;
				}
			else
				{
				rec = null;
				}
			
			if(		rec==null ||
					(ssr!=null && ssr.getSequenceIndex()!=rec.getReferenceIndex()))
				{
				if(ssr!=null)
					{
					aggregator.finish(array);
					// dump data
					int start0=(interval==null?0:interval.getStart());
					boolean header_printed=false;
					
					
					while(start0 < array.length)
						{
						if(interval!=null)
							{
							if(!interval.getContig().equals(ssr.getSequenceName())) break;//
							if(start0> interval.getEnd()) break;
							
							if(start0+window_span < interval.getStart())
								{
								start0+=win_shift;
								continue;
								}
							}
						
						if(!this.bedGraph && !header_printed)
							{
							pw.println(
			 						"fixedStep chrom="+ssr.getSequenceName()+
			 						" start="+(start0+1)+
			 						" step="+this.win_shift +" span="+ this.window_span
			 						);
							header_printed=true;
							}
						/* 
						 * http://genome.ucsc.edu/goldenPath/help/wiggle.html
						   Wiggle track data values can be integer or real, positive or negative values.
						   Chromosome positions are specified as 1-relative.
						   For a chromosome of length N, the first position is 1 and the last position is N. Only positions specified have data. Positions not specified do not have data and will not be graphed. 
						 */
		 				
						final double percentile_value = 
							percentile.evaluate(
									array,
									start0,
									Math.min(this.window_span,array.length-start0)
									);
						if(this.bedGraph)
							{
							pw.print(ssr.getSequenceName());
							pw.print('\t');
							pw.print(start0);
							pw.print('\t');
							pw.print(start0+this.window_span);
							pw.print('\t');
							}
						
						pw.printf(this.printfFormat,percentile_value);
	 					pw.print('\n');
		 				
		 				if(pw.checkError()) break;
		 				start0 += this.win_shift;
						}
					array = null;
					System.gc();
					ssr = null;
					}
				if(rec==null) break;
				if(pw.checkError()) break;
				}
			if(ssr==null)
				{
				System.gc();
				ssr=dict.getSequence(rec.getReferenceIndex());
				Objects.requireNonNull(ssr);
				LOG.info("Allocating int["+ssr.getSequenceLength()+"]");
				array=new int[ssr.getSequenceLength()];
				LOG.info("Allocating : Done.");
				Arrays.fill(array, 0);
				}
			aggregator.visit(array, rec);
			
			}
		progess.finish();
		iter.close();
		pw.flush();
		}
	
	@Override
	public int doWork(final List<String> args) {
			if(this.win_shift<=0) {
				LOG.error("window shift<=0");
				return -1;
			}
			if(this.window_span<=0) {
				LOG.error("window size<=0");
				return -1;
			}
			final Interval interval;
			PrintWriter pw = null;
			CloseableIterator<SAMRecord> samRecordIterator = null;
			final List<SamReader> samReaders = new ArrayList<>();
			final List<CloseableIterator<SAMRecord>> merginIterators= new ArrayList<>();
			try
				{
				final SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(htsjdk.samtools.ValidationStringency.LENIENT);
				
				if(args.isEmpty())
					{
					if(!StringUtil.isBlank(region_str)) {
						LOG.error("cannot specify region for stdin");
						return -1;
						}
					interval = null;
					samReaders.add( srf.open(SamInputResource.of(stdin())));
					samRecordIterator = samReaders.get(0).iterator();
					}
				else if(args.size()==1 && !args.get(0).endsWith(".list"))
					{
					samReaders.add( srf.open(SamInputResource.of(new File(args.get(0)))));
					if(StringUtil.isBlank(this.region_str))
						{
						samRecordIterator = samReaders.get(0).iterator();
						interval = null;
						}
					else
						{
						interval = new IntervalParser(samReaders.get(0).getFileHeader().getSequenceDictionary()).
							setContigNameIsWholeContig(true).
							parse(region_str);
						if(interval==null) 
							{
							LOG.error("Cannot parse interval "+this.region_str);
							return -1;
							}
						LOG.debug("interval "+interval);
						samRecordIterator = samReaders.get(0).query(
								interval.getContig(),
								interval.getStart(),
								interval.getEnd(),
								false
								);
						}
					}
				else
					{
					final List<File> samFiles;
					if(args.size()==1 && args.get(0).endsWith(".list"))
						{
						samFiles = IOUtils.unrollFile(new File(args.get(0)));
						}
					else
						{
						samFiles = args.stream().map(S->new File(S)).collect(Collectors.toList());
						}
					if(samFiles.isEmpty()) {
						LOG.error("No Input SAM file");
						return -1;
						}
					final SAMSequenceDictionary dict0 = SAMSequenceDictionaryExtractor.extractDictionary(samFiles.get(0));
					if(dict0==null) throw new JvarkitException.DictionaryMissing(samFiles.get(0).getPath());
					samFiles.stream().forEach(F->{
						final SAMSequenceDictionary dicti = SAMSequenceDictionaryExtractor.extractDictionary(F);
						if(dicti==null) throw new JvarkitException.DictionaryMissing(F.getPath());
						if(!SequenceUtil.areSequenceDictionariesEqual(dicti, dict0)) {
							throw new JvarkitException.DictionariesAreNotTheSame(dict0,dicti);
							}
						});
					for(final File bamFile: samFiles)
						{
						LOG.info("opening "+bamFile);
						samReaders.add(srf.open(bamFile));
						}
					final SamFileHeaderMerger mergedheader = new SamFileHeaderMerger(
							SAMFileHeader.SortOrder.coordinate,
							samReaders.stream().map(SR->SR.getFileHeader()).collect(Collectors.toList()),
							false
							);
					final Map<SamReader,CloseableIterator<SAMRecord>> reader2iter= new HashMap<>();
							
					if(StringUtil.isBlank(this.region_str))
						{
						for(final SamReader sr:samReaders)
							{
							reader2iter.put(sr, sr.iterator());
							}
						interval = null;
						}
					else
						{
						interval = new IntervalParser(dict0).
							setContigNameIsWholeContig(true).
							parse(region_str);
						if(interval==null) 
							{
							LOG.error("Cannot parse interval "+this.region_str);
							return -1;
							}
						LOG.info("interval :"+interval);
						for(final SamReader sr:samReaders)
							{
							reader2iter.put(sr, sr.query(
									interval.getContig(),
									interval.getStart(),
									interval.getEnd(),
									false
									));
							}
						}
					merginIterators.addAll(reader2iter.values());
					samRecordIterator = new MergingSamRecordIterator(mergedheader, reader2iter, true);
					}
				
				for(final SamReader sr:samReaders)
					{
					if(sr.getFileHeader().getSortOrder()!=SAMFileHeader.SortOrder.coordinate) {
						LOG.error("one of your bam input is not sorted on coordinate");
						return -1;
						}
					}
				pw = openFileOrStdoutAsPrintWriter(this.outputFile);
				
				run(
					pw,
					samRecordIterator,
					samReaders.get(0).getFileHeader().getSequenceDictionary(),
					interval
					);
				samRecordIterator.close();
				samRecordIterator=null;
				CloserUtil.close(samReaders);
				samReaders.clear();
				pw.flush();
				return RETURN_OK;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(merginIterators);
				CloserUtil.close(samRecordIterator);
				CloserUtil.close(samReaders);
				CloserUtil.close(pw);
				pw=null;
				}
			}
	
	public static void main(final String[] args)
		{
		new Bam2Wig().instanceMainWithExit(args);
		}

	}
