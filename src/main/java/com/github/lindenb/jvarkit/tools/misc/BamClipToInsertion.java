
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
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class BamClipToInsertion
	extends AbstractBamClipToInsertion
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AbstractBamClipToInsertion.class);
	private static class Base
		{
		CigarOperator op=null;
		int refIndex=-1;
		int readIndex=-1;
		int countInsertion=0;
		}
	
	private static class SamAndCigar
		{
		SAMRecord rec;
		List<Base> bases;
		boolean containsIorS=false;
		
		SamAndCigar(final SAMRecord rec) {
			this.rec=rec;
			this.bases = new ArrayList<>(rec.getReadLength()+1);
			Cigar cigar = rec.getCigar();
			if(cigar==null) throw new RuntimeException("no cigar in "+rec.getSAMString());
			int readpos=0;
			int refpos= rec.getUnclippedStart();
			for(int i=0;i < cigar.numCigarElements();++i)
				{
				final CigarElement ce = cigar.getCigarElement(i);
				final CigarOperator operator = ce.getOperator();
				switch(operator)
					{
					case S: case I: containsIorS = true; break;
					default:break;
					}
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
			final SAMReadGroupRecord g0 = this.rec.getReadGroup();
			if(g0==null) return ;
			final SAMReadGroupRecord g1 = other.rec.getReadGroup();
			if(g1==null || !g1.equivalent(g0)) return ;
			
			int index3 =this.bases.size();
			if(index3>1 && this.bases.get(index3-1).op.equals(CigarOperator.S))
				{
				index3--;
				}
			if(index3!=this.bases.size())
				{
				Base base3a = this.bases.get(index3);
				/* find equivalent position on 'other' */
				int i=0;
				while(i< other.bases.size())
					{
					Base base3b = other.bases.get(i);
					/* same ref */
					if( base3b.refIndex == base3a.refIndex)  break;
					++i;
					}
				/* zip over cigar while this is SOFT CLIP and other is INSERTION */
				while(i< other.bases.size() &&
					  index3< this.bases.size() &&
					  other.bases.get(i).op.equals(CigarOperator.I) &&
					  this.bases.get(index3).op.equals(CigarOperator.S)
					  )
					{
					LOG.info("Found "+this.rec.getSAMString() +"\n"+other.rec.getSAMString());
					this.bases.get(index3).countInsertion++;
					++i;
					++index3;
					}
				}
			}
		
		public SAMRecord getSAMRecord() {
			return rec;
			}
		
		public SAMRecord build()
			{
			boolean newCigar=false;
			for(int i=0;i< this.bases.size();++i)
				{
				Base base = this.bases.get(i);
				if(base.countInsertion>0)
					{
					base.op=CigarOperator.I;
					newCigar=true;
					}
				}
			if(!newCigar) return rec;
			LOG.info("NEW!");
			List<CigarElement> cigarElements = new ArrayList<>();
			int i=0;
			while(i< this.bases.size())
				{
				int j=i+1;
				while(j< this.bases.size() && 
					this.bases.get(j).op.equals(this.bases.get(i).op)	
					)
					{
					++j;
					}
				cigarElements.add(new CigarElement(j-i, this.bases.get(i).op));
				i=j;
				}
			Cigar cigar = new Cigar(cigarElements);
			rec.setCigar(cigar);
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
					//ignore unmapped reads
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
				final SamAndCigar sac = new  SamAndCigar(rec);
				if(!sac.containsIorS )
					{
					sfw.addAlignment(rec);
					continue;
					}
				buffer.add(sac);
				
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
						sfw.addAlignment(ri.build());
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
