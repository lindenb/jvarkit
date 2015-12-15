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

import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SortSamRefName extends AbstractSortSamRefName
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SortSamRefName.class);
	
	
	private static class RefNameComparator implements
		Comparator<SAMRecord>
		{
		private SAMRecordQueryNameComparator nameCmp=new SAMRecordQueryNameComparator();
		
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
	protected Collection<Throwable> call(String inputName) throws Exception {	
		SamReader in=null;
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		CloseableIterator<SAMRecord> iter2=null;
		SortingCollection<SAMRecord> sorter=null;
		try
			{
			in  = openSamReader(inputName);
			SAMFileHeader header= in.getFileHeader();
			
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(header);
			final RefNameComparator refNameComparator=new RefNameComparator();
			sorter =SortingCollection.newInstance(
					SAMRecord.class,
					bamRecordCodec,
					refNameComparator,
					super.getMaxRecordsInRam(),
					getTmpDirectories()
					);
			sorter.setDestructiveIteration(true);
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			iter = in.iterator();
			while(iter.hasNext())
				{
				sorter.add( progress.watch(iter.next()));
				}
			in.close();in=null;
			sorter.doneAdding();
			
			SAMFileHeader header2=header.clone();
			header2.addComment(getName()+" "+getVersion()+" "+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			out = super.openSAMFileWriter(header2, true);
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
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter2);
			CloserUtil.close(iter);
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		}
	
	public static void main(String[] args) throws IOException
		{
		new SortSamRefName().instanceMainWithExit(args);
		}
	}
