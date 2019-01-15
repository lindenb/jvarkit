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
* 2016 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;


/**

BEGIN_DOC

## Examples

```
$ java -jar dist/concatsam.jar dir/*.bam
$ java -jar dist/concatsam.jar --merge dir/*.bam
$ java -jar dist/concatsam.jar --region 'ref:100-200' dir/*.bam
```


END_DOC

 */
@Program(name="concatsam",
	description="concat sam files",
	keywords={"sam","bam"})
public class ConcatSam extends Launcher
	{
	private static final Logger LOG = Logger.build(ConcatSam.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","--region","--interval"},description="Limit analysis to this interval. "+IntervalParser.OPT_DESC)
	private String region_str=null;
	@Parameter(names={"-merge","--merge"},description="Don't really concatenate one sam after the other, use a htsjdk.samtools.MergingSamRecordIterator. Similar to Picard MergeSamFiles" )
	private boolean merging=false;

	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	/** a CloseableIterator<SAMRecord> with a getFileHeader function */
	public static interface ConcatSamIterator
		extends CloseableIterator<SAMRecord>
		{
		public SAMFileHeader getFileHeader();
		}
	
	/** implementation of ConcatSamIterator */
	private static class ConcatSamIteratorImpl
	extends AbstractIterator<SAMRecord>
	implements ConcatSamIterator
		{
		private final List<SamReader> samReaders = new ArrayList<>();
		private final List<CloseableIterator<SAMRecord>> merginIterators= new ArrayList<>();
		private CloseableIterator<SAMRecord> delegate=null;
		private SAMFileHeader header;
		private boolean concatenate=false;
		@Override
		protected SAMRecord advance() {
			for(;;) {
				if(this.delegate.hasNext()) return delegate.next();
				if(!this.concatenate || this.merginIterators.isEmpty())
					{
					return null;
					}
				//close current iterator
				CloserUtil.close(this.delegate);
				//pop front from the list of iterators
				this.delegate = merginIterators.remove(0);
				}			
			}
		@Override
		public SAMFileHeader getFileHeader() {
			return this.header;
			}
		
		@Override
		public void close() {
			CloserUtil.close(this.delegate);
			CloserUtil.close(this.merginIterators);
			this.merginIterators.clear();
			CloserUtil.close(this.samReaders);
			this.samReaders.clear();
			}
		
		}

	
	public static class Factory
		{
		private final List<String> intervalStrList = new ArrayList<>();
		private SamReaderFactory samReaderFactory;
		private boolean enableUnrollList=true;
		private boolean concatenate=false;
		
		public Factory() {
			this.samReaderFactory = SamReaderFactory.
					makeDefault().
					validationStringency(htsjdk.samtools.ValidationStringency.LENIENT);
			}
		
		public Factory setEnableUnrollList(boolean enableUnrollList) {
			this.enableUnrollList = enableUnrollList;
			return this;
			}
		
		public boolean isEnableUnrollList() {
			return enableUnrollList;
			}
		
		public Factory setConcatenate(boolean concatenate) {
			this.concatenate = concatenate;
			return this;
		}
		
		public boolean isConcatenate() {
			return concatenate;
			}
		
		
		
		@Deprecated
		/* deprecated use addInterval */
		public Factory setInterval(final String intervalStr) {
			return addInterval(intervalStr);
			}
		
		/** add this interval string; Will be ignored if it is blank */
		public Factory addInterval(final String intervalStr) {
			if(!StringUtil.isBlank(intervalStr)) this.intervalStrList.add(intervalStr);
			return this;
			}

		public Factory setSamReaderFactory(final SamReaderFactory samReaderFactory) {
			this.samReaderFactory = samReaderFactory;
			return this;
			}
		
		public SamReaderFactory getSamReaderFactory() {
			return samReaderFactory;
			}
		
		
		private QueryInterval[] getRegionsAsQueryIntervalArray(final SamReader r) 
			{
			final List<QueryInterval> queryIntervals=new ArrayList<>();
			for(final String intervalStr:this.intervalStrList)
				{
				if(StringUtil.isBlank(intervalStr)) {
					throw new IllegalStateException("empty interval string in "+this.intervalStrList);
				}
				final SAMSequenceDictionary dict =r.getFileHeader().getSequenceDictionary();
				if(dict==null) throw new JvarkitException.BamDictionaryMissing(r.getResourceDescription());
				final Interval i= new IntervalParser(dict).
					setContigNameIsWholeContig(true).
					parse(intervalStr);
				if(i==null) throw new IllegalArgumentException("Cannot parse interval "+intervalStr);
				final int referenceIndex = dict.getSequenceIndex(i.getContig());
				if(referenceIndex<0) throw new IllegalArgumentException("tid<0 ??! for "+i);
				queryIntervals.add(new QueryInterval(referenceIndex, i.getStart(), i.getEnd()));
				}
			Collections.sort(queryIntervals);
			return queryIntervals.toArray(QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()])));
			}

		
		public ConcatSamIterator open(final List<String> args) throws IOException {
			final ConcatSamIteratorImpl myIter=new ConcatSamIteratorImpl();
			myIter.concatenate = this.concatenate;
			final SamReaderFactory srf= getSamReaderFactory();

			
			if(args.isEmpty())
				{
				LOG.info("reading from stdin");
				if(!this.intervalStrList.isEmpty()) {
					throw new SAMException("cannot specify regions for stdin ("+this.intervalStrList+")");
					}
				final SamReader reader = srf.open(SamInputResource.of(System.in));
				myIter.samReaders.add(reader);
				myIter.delegate =  reader.iterator();
				myIter.hasNext();//for try reading something to get an error if something is wrong
				myIter.header = reader.getFileHeader();
				}
			else if(args.size()==1 && !(args.get(0).endsWith(".list") && this.isEnableUnrollList()))
				{
				final SamReader reader =  srf.open(SamInputResource.of(args.get(0)));
				myIter.samReaders.add(reader);
				myIter.header = reader.getFileHeader();
				if(this.intervalStrList.isEmpty())
					{
					myIter.delegate = reader.iterator();
					}
				else
					{
					final QueryInterval intervals[] = getRegionsAsQueryIntervalArray(reader);
					myIter.delegate = reader.query( intervals, false);
					}
				}
			else
				{
				final List<String> samFiles;
				if(args.size()==1 && args.get(0).endsWith(".list") && this.isEnableUnrollList())
					{
					samFiles = new ArrayList<>(Files.lines(new File(args.get(0)).toPath()).
							filter(S->!(S.trim().isEmpty() || S.startsWith("#"))).
							collect(Collectors.toSet())
							);
					}
				else
					{
					samFiles = args;
					}
				
				
				for(final String bamFile: samFiles)
					{
					myIter.samReaders.add(srf.open(SamInputResource.of(bamFile)));
					}
				
				if(myIter.samReaders.isEmpty()) {
					throw new SAMException("No Input SAM file");
					}
				
				final SAMSequenceDictionary dict0 = myIter.samReaders.get(0).getFileHeader().getSequenceDictionary();
				if(dict0==null) throw new JvarkitException.DictionaryMissing(samFiles.get(0));
				myIter.samReaders.stream().
					skip(1L).
					forEach(F->{
					final SAMSequenceDictionary dicti = F.getFileHeader().getSequenceDictionary();
					if(dicti==null) throw new JvarkitException.DictionaryMissing("source:"+F.getResourceDescription());
					if(!SequenceUtil.areSequenceDictionariesEqual(dicti, dict0)) {
						myIter.close();
						throw new JvarkitException.DictionariesAreNotTheSame(dict0,dicti);
						}
					});
				
				final SamFileHeaderMerger mergedheader = new SamFileHeaderMerger(
						(this.concatenate?SAMFileHeader.SortOrder.unsorted:SAMFileHeader.SortOrder.coordinate),
						myIter.samReaders.stream().map(SR->SR.getFileHeader()).collect(Collectors.toList()),
						false
						);
				
				final Map<SamReader,CloseableIterator<SAMRecord>> reader2iter= new HashMap<>();
						
				if(this.intervalStrList.isEmpty())
					{
					for(final SamReader sr:myIter.samReaders)
						{
						reader2iter.put(sr, sr.iterator());
						}
					}
				else
					{
					for(final SamReader sr:myIter.samReaders)
						{
						final QueryInterval intervals[] = getRegionsAsQueryIntervalArray(sr);
						reader2iter.put(sr, sr.query(intervals,false));
						}
					}
				myIter.merginIterators.addAll(reader2iter.values());
				myIter.header = mergedheader.getMergedHeader();
				
				if(this.concatenate) {
					myIter.delegate = myIter.merginIterators.remove(0);
					}
				else
					{
					@SuppressWarnings("resource")
					final MergingSamRecordIterator mergedIter = new MergingSamRecordIterator(mergedheader, reader2iter, true);
					myIter.delegate = mergedIter;
					}
				}
			return myIter;
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
	
		SAMFileWriter out=null;
		ConcatSamIterator iter = null;
		try
			{
			final Factory factory = new Factory().
					setConcatenate(!this.merging)
					;
			if(!StringUtil.isBlank(this.region_str)) {
				factory.addInterval(region_str);
				}
			
			iter = factory.open(args);
			
			out = this.writingBamArgs.openSAMFileWriter(
					outputFile,
					iter.getFileHeader(),
					true);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(iter.getFileHeader()).logger(LOG);
			while(iter.hasNext())
				{
				out.addAlignment(progress.watch(iter.next()));
				}
			iter.close();iter=null;
			out.close();out=null;
			progress.finish();
			return RETURN_OK;
			}
		catch(final Exception err) 
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(iter);
			}
		}
	

	public static void main(final String[] args) {
		new ConcatSam().instanceMainWithExit(args);
		}
	}
