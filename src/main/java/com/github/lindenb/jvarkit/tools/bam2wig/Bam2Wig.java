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
package com.github.lindenb.jvarkit.tools.bam2wig;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
BEGIN_DOC

## Example
the input file

```bash
java -jar dist/bam2wig.jar -w 1 -s 3 -i -t -L OFF examples/toy.bam
```

```
track type=wiggle_0 name="__REPLACE_WIG_NAME__" description="__REPLACE_WIG_DESC__"
fixedStep chrom=ref start=7 step=3 span=1
1
3
3
3
1
1
0
0
1
0
2
2
1
fixedStep chrom=ref2 start=1 step=3 span=1
1
2
3
4
5
6
6
5
4
3
3
```

END_DOC
 */
@Program(name="bam2wig",
description="Bam to fixedStep Wiggle converter. Parses the cigar String to get the depth. Memory intensive: must alloc sizeof(int)*size(chrom)",
keywords={"bam","wig","wiggle"}
)
public class Bam2Wig extends Launcher
	{
	private static final Logger LOG = Logger.build(Bam2Wig.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-t","--header"},description="print a UCSC custom track header")
	private boolean custom_track = false;

	@Parameter(names={"-s","--windowShift"},description="window shift")
	private int WINDOW_SHIFT = 25 ;

	@Parameter(names={"-w","--windowSize"},description="window size")
	private int WINDOW_SIZE = 100 ;

	@Parameter(names={"-i","--integer"},description="cast to integer")
	private boolean cast_to_integer = false;

	@Parameter(names={"-g","--zerolength"},description="minimal zero-coverage length before writing a new header")
	private int min_gap = 200 ;

	@Parameter(names={"-d","--mindepth"},description="minimal depth before setting depth to zero")
	private int min_depth = 0 ;

	@Parameter(names={"--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamFilterParser.ACCEPT_ALL;

	
	public Bam2Wig()
		{
		
		}
	

	
	
	private PrintWriter pw = null;
	
	private void run(final SamReader sfr)
		{
		SAMSequenceRecord prev_ssr=null;
		SAMSequenceDictionary dict=sfr.getFileHeader().getSequenceDictionary();
		SAMRecordIterator iter=sfr.iterator();
		int array[]=null;
		SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(dict);
		if(custom_track)
			{
			pw.println("track type=wiggle_0 name=\"__REPLACE_WIG_NAME__\" description=\"__REPLACE_WIG_DESC__\"");
			}
		
		for(;;)
			{
			SAMRecord rec=null;
			if(iter.hasNext())
				{
				rec=iter.next();
				progess.watch(rec);
				if(rec.getReadUnmappedFlag()) continue;
				if(samRecordFilter.filterOut(rec)) continue;
				}
			
			if(		rec==null ||
					(prev_ssr!=null && prev_ssr.getSequenceIndex()!=rec.getReferenceIndex()))
				{
				if(prev_ssr!=null)
					{
					int start0=0;
					while(start0< array.length && array[start0]==0)
						{
						++start0;
						}
					
					int end0 = prev_ssr.getSequenceLength();
					while(end0 >0 && array[end0-1]==0)
						{
						--end0;
						}
					
					int last_non_zero_pos0=start0;
					int num_zero_regions_skipped=0;
					boolean need_print_header=true;
					while(start0 < end0)
						{
						/* 
						 * http://genome.ucsc.edu/goldenPath/help/wiggle.html
						   Wiggle track data values can be integer or real, positive or negative values.
						   Chromosome positions are specified as 1-relative.
						   For a chromosome of length N, the first position is 1 and the last position is N. Only positions specified have data. Positions not specified do not have data and will not be graphed. 
						 */
		 				
		 					int n=0;
		 					double sum=0;
		 					for(int j=0;j< WINDOW_SIZE && start0+j< array.length;++j)
		 						{
		 						sum+=array[start0+j];
		 						n++;
		 						}
		 					
		 					if(sum/n < min_depth)
		 						{
		 						sum=0;
		 						}
		 					
		 					
		 					if(sum==0)
		 						{
		 						start0+=WINDOW_SHIFT;
		 						num_zero_regions_skipped++;
		 						continue;
		 						}
		 					else
		 						{
		 						if((start0-last_non_zero_pos0)<= min_gap)
		 							{
		 							for(int r=0;r < num_zero_regions_skipped;++r)
		 								{
		 								pw.println(0);
		 								}
		 							}
		 						else
		 							{
		 							need_print_header=(num_zero_regions_skipped>0);
		 							}
		 						
		 						last_non_zero_pos0=start0;
		 						num_zero_regions_skipped=0;
		 						}
		 					
		 					if(need_print_header)
			 					{
		 						need_print_header=false;
		 						pw.println(
				 						"fixedStep chrom="+prev_ssr.getSequenceName()+
				 						" start="+(start0+1)+
				 						" step="+WINDOW_SHIFT +" span="+ WINDOW_SIZE);
			 					}
		 					
		 					if(cast_to_integer)
		 						{
		 						pw.println((int)(sum/n));
		 						}
		 					else
		 						{
		 						pw.println((float)(sum/n));
		 						}
		 					
		 				if(pw.checkError()) break;
		 				start0+=WINDOW_SHIFT;
						}
					array=null;
					System.gc();
					prev_ssr=null;
					}
				if(rec==null) break;
				if(pw.checkError()) break;
				}
			if(prev_ssr==null)
				{
				prev_ssr=dict.getSequence(rec.getReferenceIndex());
				LOG.info("Allocating int["+prev_ssr.getSequenceLength()+"]");
				array=new int[prev_ssr.getSequenceLength()];
				LOG.info("Allocating : Done.");
				Arrays.fill(array, 0);
				}
			
			final Cigar cigar=rec.getCigar();
			if(cigar==null) continue;
    		int refpos1=rec.getAlignmentStart();
    		for(CigarElement ce:cigar.getCigarElements())
    			{
    			final CigarOperator op = ce.getOperator();
    			if(op.consumesReferenceBases())
    				{
    				if(op.consumesReadBases())
    					{
    					for(int i=0;i< ce.getLength() && refpos1<= array.length;++i)
			    			{
							if(refpos1>= 1 && refpos1<=array.length)
								{
								array[refpos1-1]++;
								}
							refpos1++;
		    				}
    					}
    				else
    					{
    					refpos1+=ce.getLength();
    					}
    				}    				
    			}
			}
		progess.finish();
		iter.close();
		pw.flush();
		}
	
	@Override
	public int doWork(final List<String> args) {
			SamReader in=null;
			try
				{
				this.pw = openFileOrStdoutAsPrintWriter(this.outputFile);
				in = openSamReader(oneFileOrNull(args));
				run(in);
				pw.flush();
				return RETURN_OK;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(in);
				CloserUtil.close(pw);
				pw=null;
				in=null;
				}
			}
		
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Bam2Wig().instanceMainWithExit(args);
		}

	}
