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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;


/**
 * 
 * Biostar234230
 *
 */
public class Biostar234230 extends AbstractBiostar234230
	{
	private static class SlidingWindow
		{
		final int start;
		final int end;
		int pairs_in_window=0;
		int pairs_over_window=0;
		int pairs_partial_overlap=0;
		
		SlidingWindow(int start,int end) {
		this.start = start;
		this.end=end;
		}
		
		void print(PrintWriter out ,final String contig) {
			out.print(contig);
			out.print("\t");
			out.print( start);
			out.print("\t");
			out.print( end);
			out.print("\t");
			out.print( pairs_in_window);
			out.print("\t");
			out.print( pairs_over_window);
			out.print("\t");
			out.print( pairs_partial_overlap);
			out.println();
			}
		
		void watch(int chromStart,int chromEnd) {
			if(chromStart>=start && chromEnd<=end) {
				++pairs_in_window;
				}
			else if(chromEnd<this.start || chromStart>this.end)
				{
				//nothing
				}
			else if(chromStart<this.start && chromEnd > this.end)
				{
				++pairs_over_window;
				}
			else
				{
				++pairs_partial_overlap;
				}
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		
		
		if(super.windowSize<=0)
			{
			return wrapException("Bad window size.");
			}
		if(super.windowShift<=0)
			{
			return wrapException("Bad window shift.");
			}
		
		
		SamReader in = null;
		PrintWriter out=null;
		SAMRecordIterator iter=null;
		final List<SlidingWindow> buffer= new ArrayList<>();
		try
			{
			
			in = super.openSamReader(inputName);
			final SAMFileHeader header = in.getFileHeader();
			if(header.getSortOrder()!=SAMFileHeader.SortOrder.coordinate) {
			return wrapException("Bam is not sorted on coordinate");
			}
			out = super.openFileOrStdoutAsPrintWriter();
			
			out.print("#contig");
			out.print("\t");
			out.print("start");
			out.print("\t");
			out.print("end");
			out.print("\t");
			out.print("pairs_in_window");
			out.print("\t");
			out.print("pairs_over_window");
			out.print("\t");
			out.print("pairs_partial_overlap");
			out.println();
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			iter = in.iterator();
			String prev_contig=null;
			for(;;)
				{
				final SAMRecord rec = (iter.hasNext()?progress.watch(iter.next()):null);
				if(rec!=null )
					{
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(rec.getMateAlignmentStart()> rec.getAlignmentStart()) continue;
					}
				if(rec==null || !rec.getContig().equals(prev_contig))
					{
					for(final SlidingWindow w:buffer)
						{
						w.print(out,prev_contig);
						}
					buffer.clear();
					if(rec==null) break;
					prev_contig = rec.getContig();
					}
				final int fragStart = rec.getMateAlignmentStart();
				final int fragEnd = rec.getAlignmentEnd();
				while(!buffer.isEmpty() && buffer.get(0).end < fragStart)
					{
					buffer.remove(0).print(out,rec.getContig());;
					}
				final int winStart = fragStart - fragStart%super.windowSize;
				final int winEnd = fragEnd - fragEnd%super.windowSize;
				
				for(int winPos = winStart;winPos<=winEnd;winPos+=super.windowShift) {
					if(buffer.isEmpty() || buffer.get(buffer.size()-1).start < winPos)
						{
						buffer.add(new SlidingWindow(winPos,winPos+super.windowSize));
						}
					}
				for(SlidingWindow w:buffer)
					{
					w.watch(fragStart,fragEnd);
					}				
				}
			
				
			progress.finish();
			out.flush();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(in);
			CloserUtil.close(out);
			}	
		}
	
	public static void main(String[] args) {
		new Biostar234230().instanceMainWithExit(args);
		}
	}
