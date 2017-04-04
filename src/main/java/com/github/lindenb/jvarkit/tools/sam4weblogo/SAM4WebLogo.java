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
package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

@Program(name="sam4weblogo",
	description="Sequence logo for different alleles or generated from SAM/BAM ",
	biostars=73021,
	keywords={"sam","bam","visualization"}
	)
public class SAM4WebLogo extends Launcher
	{
	private static final Logger LOG = Logger.build(SAM4WebLogo.class).make();
	
	@Parameter(names={"-c","--clipped"},description="Use Soft Clipped Bases")
	private boolean useClip = false;

	@Parameter(names={"-r","--region"},description="Region to observe: chrom:start-end",required=true)
	private String regionStr = null;

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@Override
	public int doWork(final List<String> args) {
		
    	if(this.regionStr==null)
			{
			throw new JvarkitException.CommandLineError("Undefined interval.");
			}
    	final String inputName=oneFileOrNull(args);
    
		final Interval interval=parseInterval(this.regionStr);
		if(interval==null)
			{
			throw new JvarkitException.CommandLineError("Bad interval "+this.regionStr);
			}
		
		PrintWriter out=null;
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		boolean warningInsertion=false;
		try {
			samReader = openSamReader(inputName);
			out =  openFileOrStdoutAsPrintWriter(this.outputFile);
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
                if(useClip) {
                    if(rec.getUnclippedEnd() < interval.getStart() ) continue;
                    if(rec.getUnclippedStart() > interval.getEnd() ) continue;
                } else
	                {
	                if(rec.getAlignmentEnd() < interval.getStart() ) continue;
	                if(rec.getAlignmentStart() > interval.getEnd() ) continue;
	                }
                Cigar cigar=rec.getCigar();
                if(cigar==null) continue;
                byte bases[]=rec.getReadBases();
                
                StringBuilder seq=new StringBuilder(interval.length());
                int readPos=0;
                int refPos=rec.getUnclippedStart();
                for(int i=0;i< cigar.numCigarElements();++i)
                	{
                	final CigarElement ce=cigar.getCigarElement(i);
                	final CigarOperator op=ce.getOperator();
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
	       out.flush();
	       out.close();out=null;
	       
	       if(warningInsertion)
	        	{
	        	LOG.warn("Some reads contained insertions.");
	        	}
	        return 0;
			} 
		catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
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
