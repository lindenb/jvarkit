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




*/
package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.PrintWriter;
import java.util.Collection;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAM4WebLogo extends AbstractSAM4WebLogo
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SAM4WebLogo.class);

    @Override
    protected Collection<Throwable> call(String inputName) throws Exception {
    
    	if(super.regionStr==null)
			{
			return wrapException("Undefined interval.");
			}
    	
    
		final Interval interval=parseInterval(this.regionStr);
		if(interval==null)
			{
			return wrapException("Bad interval "+this.regionStr);
			}
		
		PrintWriter out= openFileOrStdoutAsPrintWriter();
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		boolean warningInsertion=false;
		try {
			samReader = openSamReader(inputName);
		
			if(samReader.hasIndex())
					{
					iter=samReader.queryOverlapping(
							interval.getContig(),
							interval.getStart(),
							interval.getEnd()
							);
					}
			else
					{
					iter=samReader.iterator();
					}
			
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
	       while(iter.hasNext())
                {
                SAMRecord rec=iter.next();
                progress.watch(rec);
                if(rec.getReadUnmappedFlag()) continue;
                if(!rec.getReferenceName().equals(interval.getContig())) continue;
                if(rec.getAlignmentEnd() < interval.getStart() ) continue;
                if(rec.getAlignmentStart() > interval.getEnd() ) continue;
                Cigar cigar=rec.getCigar();
                if(cigar==null) continue;
                byte bases[]=rec.getReadBases();
                
                StringBuilder seq=new StringBuilder(interval.length());
                int readPos=0;
                int refPos=rec.getUnclippedStart();
                for(int i=0;i< cigar.numCigarElements();++i)
                	{
                	CigarElement ce=cigar.getCigarElement(i);
            		CigarOperator op=ce.getOperator();
                	switch(op)
                		{
                		case P:break;
                		case I:
                			{
                			warningInsertion=true;
                			break;
                			}
                		case D: case N:
                			{
		        			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
		        				{
		        				if(refPos>= interval.getStart())
		        					{
		        					seq.append('-');
		        					}
		        				refPos++;
		        				}
		        			break;
                			}
                		case H:
		        			{
		        			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
		        				{
		        				if(refPos>= interval.getStart() && useClip)
		        					{
		        					seq.append('N');
		        					}
		        				refPos++;
		        				}
		        			break;
		        			}
                		case S:
                			{
                			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
                				{
                				if(refPos>= interval.getStart() && useClip)
                					{
                					seq.append((char)bases[readPos]);
                					}
                				readPos++;
                				refPos++;
                				}
                			break;
                			}
                		case M:case X: case EQ:
                			{
                			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd() ;++j)
                				{
                				if(refPos>= interval.getStart())
                					{
                					seq.append((char)bases[readPos]);
                					}
                				readPos++;
                				refPos++;
                				}
                			break;
                			}
                		default:throw new IllegalStateException("Not handled. op:"+op);
                		}	
                	}
                if(seq.length()==0) continue;
                for(int i= interval.getStart();
                	i< (useClip?rec.getUnclippedStart():rec.getAlignmentStart());
                	++i)
                	{
                	seq.insert(0, '-');
                	}
                while(seq.length()< interval.length())
	            	{
	            	seq.append('-');
	            	}
            	out.print(">"+rec.getReadName());
            	if(rec.getReadPairedFlag())
	            	{
	            	if(rec.getFirstOfPairFlag()) out.print("/1");
	            	if(rec.getSecondOfPairFlag()) out.print("/2");
	            	}
            	out.println();
            	out.println(seq);
                }
	       progress.finish();
	       if(warningInsertion)
	        	{
	        	LOG.warn("Some reads contained insertions.");
	        	}
	        return RETURN_OK;
			} 
		catch (Exception e) {
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			out.flush();
			}
		}

private static Interval parseInterval(String reg)
	{
	try
			{	 
				int colon = reg.indexOf(':');
				if(colon<1) throw new IllegalArgumentException("bad region "+reg);
				int hyphen = reg.indexOf('-');
		
				String s=reg.substring(0,colon);
				int start= Integer.parseInt(reg.substring(colon+1,hyphen));
				int end=Integer.parseInt(reg.substring(hyphen+1));
				return new Interval(s, start, end);
			}
			catch(Exception err)
			{
				System.err.println("bad interval "+reg);
				return null;
			}
		}
public static void main(String[] args)
	{
	new SAM4WebLogo().instanceMainWithExit(args);
	}
}
