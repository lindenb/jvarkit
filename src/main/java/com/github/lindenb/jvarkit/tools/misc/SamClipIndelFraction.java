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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
 * 
BEGIN_DOC

## Example

```bash
$ samtools view -h -F3844 my.bam  | java -jar dist/samclipindelfraction.jar 

##UNMAPPED_READS=0
##MAPPED_READS=3028359
##CLIPPED_READS=1182730
##CLIPPED_READS_5_PRIME=597757
##CLIPPED_READS_3_PRIME=617399
##UNCLIPPED_READS=1845629
##COUNT_BASES=338644685
#CLIP	COUNT	FRACTION_OF_MAPPED_READS
0	1845629	0.5
1	7	1.8963724562195327E-6
2	6756	0.0018302703306027376
3	695	1.8828269386751074E-4
4	794	2.1510281860547272E-4
5	819	2.2187557737768533E-4
6	471	1.275987752684857E-4
7	447	1.210969268471616E-4
(...)
```

plotting:
```bash
$ java -jar dist/samclipindelfraction.jar |\
   grep -v "##" | cut -f1,2 | tr -d '#' > output.txt
```

then, in R:
```R
T<-read.table('output.txt',header=TRUE)
plot(T[T$CLIP>0,])
```

END_DOC
 */
@Program(
		name="Samclipindelfraction",
		description="Extract clipping/indel fraction from BAM",
		deprecatedMsg="This tool can be replace with Bioalcidaejdk",
		keywords={"sam","bam","clip"}
		)
public class SamClipIndelFraction extends Launcher
	{
	private static final Logger LOG = Logger.build(SamClipIndelFraction.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private enum Type {
		leftclip,
		rightclip,
		allclip,
		insert,
		deletion
		};
	@Parameter(names="-t",description="type")
	private Type type=Type.allclip;
		
	@Override
	public int doWork(final List<String> args)
		{
		SamReader sfr=null;
		SAMRecordIterator iter=null;
		PrintWriter pw=null;
		try
			{
			sfr = openSamReader(oneFileOrNull(args));
			pw =  super.openFileOrStdoutAsPrintWriter(outputFile);
			long total_bases_count=0L;
			long count_clipped_reads=0L;
			long count_clipped_left_reads=0L;
			long count_clipped_right_reads=0L;
			long count_unclipped_reads=0L;
			long count_unmapped_reads=0L;
			SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(sfr.getFileHeader()).logger(LOG);
			Counter<Integer> counter=new Counter<>();
			iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord record=progress.watch(iter.next());
				if(record.getReadUnmappedFlag())
					{
					++count_unmapped_reads;
					continue;
					}
				final Cigar cigar=record.getCigar();
				int left_clip_length=0;
				int right_clip_length=0;
				int deletion_N_length=0;
				int deletion_D_length=0;
				int insertion_length=0;
				
				boolean onLeft=true;
				
				for(int i=0;i<  cigar.numCigarElements();++i)
					{
					final CigarElement ce= cigar.getCigarElement(i);
					final CigarOperator op = ce.getOperator();
					
					switch(op)
						{
						case N:
							{
							onLeft=false;
							deletion_D_length+=ce.getLength();
							total_bases_count+=ce.getLength();
							break;
							}
						case D:
							{
							onLeft=false;
							deletion_N_length+=ce.getLength();
							total_bases_count+=ce.getLength();
							break;
							}
						case I:
							{
							onLeft=false;
							insertion_length+=ce.getLength();
							total_bases_count+=ce.getLength();
							break;
							}
						case S:
						case H:
							{
							if(onLeft)
								{
								if(record.getReadNegativeStrandFlag())
									{
									right_clip_length += ce.getLength();
									}
								else
									{
									left_clip_length += ce.getLength();
									}
								}
							else
								{
								if(record.getReadNegativeStrandFlag())
									{
									left_clip_length += ce.getLength();
									}
								else
									{
									right_clip_length += ce.getLength();
									}
								}
							total_bases_count+=ce.getLength();
							break;
							}
						default:
							{	
							onLeft=false;
							if(op.consumesReadBases())
								{
								total_bases_count+=ce.getLength();
								}
							break;
							}
						}
					}
				
				if( left_clip_length + right_clip_length == 0)
					{
					count_unclipped_reads++;
					}
				else
					{
					if(left_clip_length>0) count_clipped_left_reads++;
					if(right_clip_length>0) count_clipped_right_reads++;
					count_clipped_reads++;
					}
				
				switch(type)
					{
					case leftclip:  counter.incr(left_clip_length); break;
					case rightclip:  counter.incr(right_clip_length); break;
					case allclip:  counter.incr(left_clip_length+right_clip_length); break;
					case deletion: counter.incr(deletion_D_length+ deletion_N_length); break;
					case insert: counter.incr(insertion_length); break;
					default: LOG.error("Bad type: "+type);return -1;
					}
				}
			progress.finish();
			
			pw.println("##UNMAPPED_READS="+count_unmapped_reads);
			pw.println("##MAPPED_READS="+(count_clipped_reads+count_unclipped_reads));
			pw.println("##CLIPPED_READS="+count_clipped_reads);
			pw.println("##CLIPPED_READS_5_PRIME="+count_clipped_left_reads);
			pw.println("##CLIPPED_READS_3_PRIME="+count_clipped_right_reads);
			pw.println("##UNCLIPPED_READS="+count_unclipped_reads);
			pw.println("##COUNT_BASES="+total_bases_count);
			pw.print("#");
			switch(type)
				{
				case leftclip:  pw.print("CLIP_5_PRIME"); break;
				case rightclip:  pw.print("CLIP_3_PRIME"); break;
				case allclip:  pw.print("CLIP"); break;
				case deletion:  pw.print("DELETION"); break;
				case insert:  pw.print("INSERTION"); break;
				default: LOG.error("Bad type: "+type);return -1;
				}
			pw.println("\tCOUNT\tFRACTION_OF_MAPPED_READS");
			
			for(final Integer size: new TreeSet<Integer>(counter.keySet()))
				{
				pw.print(size);
				pw.print('\t');
				pw.print(counter.count(size));
				pw.print('\t');
				pw.println(counter.count(size)/(double)(count_unclipped_reads+count_unclipped_reads));
				}
			pw.flush();
			pw.close();
			pw=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(pw);
			}
		}
	
	public static void main(final String[] args) {
		new SamClipIndelFraction().instanceMainWithExit(args);

	}

}
