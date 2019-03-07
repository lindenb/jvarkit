/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
 
/**
BEGIN_DOC

## Motivation

Because gatk is buggy: http://gatkforums.broadinstitute.org/gatk/discussion/8342/duplicate-columns-in-readlengthdistribution#latest


## Input

input is a set of bam file or a file with suffix '.list' containing the path to the bam


## Example

```
$ java -jar dist/samreadlengthdistribution.jar src/test/resources/S*.bam

#ReadLength     S1      S2      S3      S4      S5
[-Inf/0[	0	0	0	0	0
[0/10[	0	0	0	0	0
[10/20[	0	0	0	0	0
[20/30[	0	0	0	0	0
[30/40[	0	0	0	0	0
[40/50[	0	0	0	0	0
[50/100[	1998	1998	1998	1998	1998
[100/150[	0	0	0	0	0
[150/200[	0	0	0	0	0
[200/250[	0	0	0	0	0
[250/300[	0	0	0	0	0
[300/350[	0	0	0	0	0
[350/400[	0	0	0	0	0
[400/450[	0	0	0	0	0
[450/500[	0	0	0	0	0
[500/1000[	0	0	0	0	0
[1000/Inf[	0	0	0	0	0
```


END_DOC
 */
@Program(name="samreadlengthdistribution",
description="Sam read length distribution",
keywords={"sam","bam","histogram"},
modificationDate="20190307"
)
public class SamReadLengthDistribution extends Launcher
	{
	private static final Logger LOG = Logger.build(SamReadLengthDistribution.class).make();
	
	private enum Method {
		SEQ_LENGTH,
		CIGAR_REF_LENGTH,
		CIGAR_PADDED_REF_LENGTH,
	}
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition=SAMRecordPartition.sample;
	@Parameter(names={"-m"},description="method, how should I get the read length ?")
	private Method method = Method.SEQ_LENGTH;
	@Parameter(names={"-w","--windows"},description=RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers ranges = new RangeOfIntegers(0,10,20,30,40,50,100,150,200,250,300,350,400,450,500,1000);

	private final Map<String,Counter<RangeOfIntegers.Range>> lengths=new TreeMap<>();

	
    public SamReadLengthDistribution()
    	{
    	}		

  
    
    private void scan(final SamReader in,String fname) throws IOException
    	{
    	final String defName = fname+"#"+this.partition.name();
    	in.getFileHeader().getReadGroups().
    		stream().
    		map(RG->this.partition.apply(RG)).
    		map(S->StringUtils.isBlank(S)?defName:S).
    		filter(S->!this.lengths.containsKey(S)).
    		forEach(S->this.lengths.put(S,new Counter<>()));
    	
    	
    	final SAMRecordIterator iter=in.iterator();
    	while(iter.hasNext())
			{
    		final SAMRecord rec = iter.next();
    		final String sampleName= this.partition.getPartion(rec,defName);
    		
    		Counter<RangeOfIntegers.Range> counter = this.lengths.get(sampleName);
    		if(counter==null) {
    			counter =  new Counter<>();
    			this.lengths.put(sampleName,counter);
    			}
    		
    		final int len;
    		switch(this.method) {
    			case SEQ_LENGTH: len = rec.getReadLength(); break;
    			case CIGAR_REF_LENGTH: {
    				if(rec.getReadUnmappedFlag()) continue;
    				final Cigar c=rec.getCigar();
    				if(c==null) continue;
    				len = c.getReferenceLength();
    				break;
    				}
    			case CIGAR_PADDED_REF_LENGTH: {
    				if(rec.getReadUnmappedFlag()) continue;
    				final Cigar c=rec.getCigar();
    				if(c==null) continue;
    				len = c.getPaddedReferenceLength();
    				break;
    				}
    			default: throw new IllegalStateException("unsupported method " +this.method);
    			}
    		
    		
    		
    		counter.incr(this.ranges.getRange(len));    		
			}
    	iter.close();
    	}
    

	@Override
	public int doWork(final List<String> inputs)
		{
		this.lengths.clear();
		SamReader r= null;
		PrintWriter out = null;
		try
			{
			final Set<String> args = IOUtils.unrollFiles(inputs);

			
			if(args.isEmpty())
				{
				r = super.openSamReader(null);
				scan(r,"<STDIN>");
				r.close();
				}
			else
				{
				for(final String filename: args)
					{
					r = super.openSamReader(filename);
					scan(r,filename);
					r.close();
					}
				}
			out=super.openFileOrStdoutAsPrintWriter(this.outputFile);
			out.print("#ReadLength");
			for(final String sample:this.lengths.keySet())
				{
				out.print("\t");
				out.print(sample);
				}
			out.println();
			
					
			
			for(final RangeOfIntegers.Range L:this.ranges.getRanges()) {
				out.print(L);
				for(final String sample:this.lengths.keySet())
					{
					final Counter<RangeOfIntegers.Range> c=this.lengths.get(sample);
					out.print("\t");
					out.print(c.count(L));
					}
				out.println();
				}
			
			
			out.flush();
			out.close();
			out = null;
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new SamReadLengthDistribution().instanceMainWithExit(args);

	}

}
