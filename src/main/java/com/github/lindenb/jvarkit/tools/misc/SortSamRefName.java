/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="sortsamrefname",description="Sort a BAM of contig and then on name")
public class SortSamRefName extends Launcher
	{
	private static final Logger LOG = Logger.build(SortSamRefName.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();


	
	
	private static class RefNameComparator implements
		Comparator<SAMRecord>
		{
		private final SAMRecordQueryNameComparator nameCmp=new SAMRecordQueryNameComparator();
		
		@Override
		public int compare(final SAMRecord o1, final SAMRecord o2)
			{
	        final int refIndex1 = o1.getReferenceIndex();
	        final int refIndex2 = o2.getReferenceIndex();
	        final int cmp = refIndex1 - refIndex2;
	        
	        if (cmp != 0)
	        	{
	        	if (refIndex1 == -1 ) return 1;
	        	if (refIndex2 == -1 ) return -1;
	            return cmp;
	        	}
			return nameCmp.compare(o1, o2);
			}
		
		}
	
	public SortSamRefName()
		{
		}

	@Override
	public int doWork(List<String> args) {
		SamReader in=null;
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		CloseableIterator<SAMRecord> iter2=null;
		SortingCollection<SAMRecord> sorter=null;
		try
			{
			in  = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header= in.getFileHeader();
			
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(header);
			final RefNameComparator refNameComparator=new RefNameComparator();
			sorter =SortingCollection.newInstance(
					SAMRecord.class,
					bamRecordCodec,
					refNameComparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpDirectories()
					);
			sorter.setDestructiveIteration(true);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			iter = in.iterator();
			while(iter.hasNext())
				{
				sorter.add( progress.watch(iter.next()));
				}
			in.close();in=null;
			sorter.doneAdding();
			
			final SAMFileHeader header2=header.clone();
			header2.addComment(getProgramName()+" "+getVersion()+" "+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			out = this.writingBamArgs.openSAMFileWriter(outputFile,header2, true);
			iter2 = sorter.iterator();
			while(iter2.hasNext())
				{
				out.addAlignment(iter2.next());
				}
			out.close();out=null;
			sorter.cleanup();
			progress.finish();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter2);
			CloserUtil.close(iter);
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		}
	
	public static void main(final String[] args) throws IOException
		{
		new SortSamRefName().instanceMainWithExit(args);
		}
	}
