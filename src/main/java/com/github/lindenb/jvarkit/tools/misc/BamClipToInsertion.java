
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
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class BamClipToInsertion
	extends AbstractBamClipToInsertion
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BamClipToInsertion.class);
	private static class Base
		{
		CigarOperator op=null;
		int refIndex=-1;
		int readIndex=-1;
		}
	
	private static class SamAndCigar
		{
		SAMRecord rec;
		List<Base> bases;
		SamAndCigar(final SAMRecord rec) {
			this.rec=rec;
			this.bases = new ArrayList<>(rec.getReadLength()+1);
			Cigar cigar = rec.getCigar();
			int readpos=0;
			int refpos= rec.getUnclippedStart();
			for(int i=0;i < cigar.numCigarElements();++i)
				{
				final CigarElement ce = cigar.getCigarElement(i);
				final CigarOperator operator = ce.getOperator();
				for(int j=0;j< ce.getLength();++j)
					{
					Base base = new Base();
					base.op = operator;
					if(base.op.consumesReadBases())
						{
						base.readIndex = readpos;
						++readpos;
						}
					else
						{
						base.readIndex = -1;
						}
					if(base.op.consumesReferenceBases()  ||
						base.op==CigarOperator.S ) //because getUnclippedStart used
						{
						base.refIndex = refpos;
						++refpos;
						}
					else
						{
						base.refIndex = -1;
						}
					this.bases.add(base);
					}
				}
			}
		void merge(final SamAndCigar other)
			{
			if(this.rec.getAlignmentEnd()< this.rec.getUnclippedEnd())
				{
				int startClip = this.bases.size();
				//scan 3' to get non clipped base
				while((startClip-1)>=0 && this.bases.get(startClip-1).op==CigarOperator.S)
					{
					--startClip;
					}
				
				}
			}
		
		public SAMRecord getSAMRecord() {
			return rec;
			}
		}
	
	public BamClipToInsertion()
		{
		}

	@Override
	public Collection<Throwable> call(final String inputName) throws Exception 
		{					SAMRecordIterator iter=null;
		SamReader sfr=null;
		SAMFileWriter sfw =null;
		try
			{			
			sfr = openSamReader(inputName);
			
			final SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				return wrapException("File header missing");
				}
			
			if(header1.getSortOrder()!=SortOrder.coordinate)
				{
				return wrapException("Input is not sorted on coordinate.");
				}
			
			final SAMFileHeader header2=header1.clone();
			header2.addComment(getName()+":"+getVersion()+":"+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			
			sfw =  openSAMFileWriter(header2, true);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1);
			iter=sfr.iterator();
			String curContig=null;
			LinkedList<SamAndCigar> buffer=new LinkedList<>();
			for(;;)
				{
				SAMRecord rec=null;
				if( iter.hasNext())
					{
					rec= progress.watch(iter.next());
					if(rec.getReadUnmappedFlag())
						{
						sfw.addAlignment(rec);
						continue;
						}
					}
				
				if(rec==null || (curContig!=null && !curContig.equals(rec.getReferenceName())))
					{
					for(final SamAndCigar r: buffer) sfw.addAlignment(r.getSAMRecord());
					buffer.clear();
					// we're done
					if(rec==null) break;
 					curContig = rec.getReferenceName();
					}
				
				buffer.add(new SamAndCigar(rec));
				
				int i=0;
				while( i < buffer.size())
					{
					final SamAndCigar ri = buffer.get(i);
					if(ri.getSAMRecord().getUnclippedEnd() < rec.getUnclippedStart())
						{
						for(int j=0;j< buffer.size();++j)
							{
							if(i==j) continue;
							ri.merge(buffer.get(j));
							}
						sfw.addAlignment(ri.getSAMRecord());
						buffer.remove(i);
						}
					else
						{
						++i;
						}
					}
				
				}
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
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
				
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamClipToInsertion().instanceMainWithExit(args);
		}
	}
