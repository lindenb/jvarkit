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

import java.io.PrintStream;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;


public class Biostar78285 extends AbstractBiostar78285
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar76892.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	
	static private class MyCommand extends AbstractBiostar78285.AbstractBiostar78285Command
		{    
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			SAMRecordIterator iter=null;
			SamReader samFileReader=null;
			PrintStream out=null;
	    	try
		    	{
	    		if(getOutputFile()==null)
	    			{
	    			out= stdout();
	    			}
	    		else
	    			{
	    			out = new PrintStream(getOutputFile());
	    			}
	    		
	    		samFileReader = openSamReader(inputName);
		    	final SAMFileHeader header=samFileReader.getFileHeader();
		    	if(header.getSortOrder()!=SortOrder.coordinate)
		    		{
		    		return wrapException("Sam file is not sorted on coordinate :"+header.getSortOrder());
		    		}
		    	SAMSequenceDictionary dict=header.getSequenceDictionary();
		    	if(dict==null)
		    		{
		    		return wrapException("SamFile doesn't contain a SAMSequenceDictionary.");
		    		}
		    	/* flag, do we saw all chromosomes in dictionary ? */
		    	boolean seen_tid[]=new boolean[dict.getSequences().size()];
		    	Arrays.fill(seen_tid, false);
		    	
		    	
	    		BitSet mapped=null;
	    		SAMSequenceRecord ssr=null;
	    		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
	    		iter=samFileReader.iterator();
		    	while(iter.hasNext())
		    		{
		    		SAMRecord rec=progress.watch(iter.next());
		    		if(rec.getReadUnmappedFlag()) continue;
		    		if(rec.isSecondaryOrSupplementary()) continue;
		    		if(rec.getDuplicateReadFlag()) continue;
		    		if(rec.getReadFailsVendorQualityCheckFlag()) continue;
		    		if(rec.getMappingQuality()==0) continue;
		    		Cigar cigar=rec.getCigar();
		    		if(cigar==null) continue;
		    		if(ssr==null || ssr.getSequenceIndex()!=rec.getReferenceIndex())
		    			{
		    			if(ssr!=null && mapped!=null)
		    				{
			    			dump(out,ssr,mapped);
			    			}
		    			ssr=dict.getSequence(rec.getReferenceIndex());
		    			if(ssr==null)
		    				{
		    				LOG.error("Sequence not in dict :"+rec);
		    				}
		    			LOG.info("allocating bitset for "+ssr.getSequenceName()+" LENGTH="+ssr.getSequenceLength());
		    			mapped=new BitSet(ssr.getSequenceLength());
		    			seen_tid[rec.getReferenceIndex()]=true;
		    			}
		    		int refpos0=rec.getAlignmentStart()-1;
		    		for(CigarElement ce:cigar.getCigarElements())
		    			{
		    			CigarOperator op=ce.getOperator();
		    			if(op.consumesReferenceBases())
		    				{	
		    				if(op.consumesReadBases())
		    					{
		    					for(int i=0;i< ce.getLength() && refpos0 +i < ssr.getSequenceLength();++i)
		    		    			{
		    						mapped.set(refpos0+i,true);
	    		    				}
		    					}
		    				refpos0 += ce.getLength();
		    				}		    				
		    			}
		    		}
		    	
		    	
		    	if(ssr!=null && mapped!=null)
					{
	    			dump(out,ssr,mapped);
	    			}
		    		
					
					
		    	/* unseen chromosomes */
		    	for(int i=0;i< seen_tid.length;++i)
		    		{
		    		if(seen_tid[i]) continue;
		    		ssr=dict.getSequence(i);
	    			out.println(ssr.getSequenceName()+"\t0\t"+ssr.getSequenceLength());
		    		}
		    	
		    	progress.finish();
		    	out.flush();
		    	out.close();
		    	return Collections.emptyList();
		    	}
	    	catch(Exception err)
	    		{
	    		return wrapException(err);
	    		}
	    	finally
	    		{
	    		CloserUtil.close(iter);
	    		CloserUtil.close(samFileReader);
	    		CloserUtil.close(out);
	    		}
	    	
	    	}
	    private void dump(PrintStream out,SAMSequenceRecord ssr,BitSet mapped)
	    	{
	    	int i=0;
	    	while(i<ssr.getSequenceLength())
	    		{
	    		if(mapped.get(i))
	    			{
	    			++i;
	    			continue;
	    			}
	    		int j=i+1;
	    		while(j<ssr.getSequenceLength() && !mapped.get(j))
	    			{
	    			++j;
	        		}
	    		out.println(ssr.getSequenceName()+"\t"+i+"\t"+j);
	    		i=j;
	    		}
	    	}
	   
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar78285().instanceMainWithExit(args);

	}

}
