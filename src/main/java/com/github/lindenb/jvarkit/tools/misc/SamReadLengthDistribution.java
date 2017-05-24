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

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/*
BEGIN_DOC

Because gatk is buggy: http://gatkforums.broadinstitute.org/gatk/discussion/8342/duplicate-columns-in-readlengthdistribution#latest

END_DOC

 */
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
 
/**
BEGIN_DOC



END_DOC
 */
@Program(name="samreadlengthdistribution",
description="Sam read length distribution",
keywords={"sam","bam","histogram"}
)
public class SamReadLengthDistribution extends Launcher
	{
	private static final Logger LOG = Logger.build(SamReadLengthDistribution.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by")
	private SAMRecordPartition partition=SAMRecordPartition.sample;
	
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
    		String sampleName= this.partition.getPartion(rec);
    		if(sampleName==null || sampleName.isEmpty()) sampleName="__UNDEFINED_"+this.partition.name().toUpperCase()+"_";
    		
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
	public int doWork(final List<String> inputs)
		{
		this.lengths.clear();
		this.max_length=0;
		SamReader r= null;
		try
			{
			final Set<String> args = IOUtils.unrollFiles(inputs);

			
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
			final PrintWriter out=super.openFileOrStdoutAsPrintWriter(this.outputFile);
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
			LOG.error(err);
			return -1;
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
