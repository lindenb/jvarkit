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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;

public class SamReadLengthDistribution extends AbstractSamReadLengthDistribution
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamReadLengthDistribution.class);

	
	private final Map<String,Counter<Integer>> lengths=new TreeMap<>();
	private int max_length=0;
	
    private SamReadLengthDistribution()
    	{
    	}		
   
    

    

    private void scan(final SamReader in) throws IOException
    	{
    	SAMRecordIterator iter=in.iterator();
    	while(iter.hasNext())
			{
    		final SAMRecord rec = iter.next();
    		final SAMReadGroupRecord srg = rec.getReadGroup();
    		String sampleName= null;
    		if(srg!=null) sampleName = srg.getSample();
    		if(sampleName==null || sampleName.isEmpty()) sampleName="__UNDEFINED_SAMPLE_";
    		
    		Counter<Integer> counter = this.lengths.get(sampleName);
    		if(counter==null) {
    			counter = new Counter<>();
    			 this.lengths.put(sampleName,counter);
    		}
    		final int len = rec.getReadLength();
    		counter.incr(len);
    		this.max_length=Math.max(max_length, len);
    		
			}
    	iter.close();
    	}
    

	
	@Override
	public Collection<Throwable> initializeKnime() {
		this.lengths.clear();
		this.max_length=0;
		return super.initializeKnime();
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		SamReader r= null;
		try
			{
			final Set<String> args = IOUtils.unrollFiles(super.getInputFiles());

			
			if(args.isEmpty())
				{
				r = super.openSamReader(null);
				scan(r);
				r.close();
				}
			else
				{
				for(final String filename: args)
					{
					r = super.openSamReader(filename);
					scan(r);
					r.close();
					}
				}
			final PrintWriter out=super.openFileOrStdoutAsPrintWriter();
			out.print("#ReadLength");
			for(final String sample:this.lengths.keySet())
				{
				out.print("\t");
				out.print(sample);
				}
			out.println();
			for(int L=0;L<=this.max_length;++L) {
				out.print(L);
				for(final String sample:this.lengths.keySet())
					{
					final Counter<Integer> c=this.lengths.get(sample);
					out.print("\t");
					out.print(c.count(L));
					}
				out.println();
				}
			
			
			out.flush();
			out.close();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{

			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamReadLengthDistribution().instanceMainWithExit(args);

	}

}
