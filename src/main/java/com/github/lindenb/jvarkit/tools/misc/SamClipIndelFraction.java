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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.util.Arrays;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
 * SamClippingFraction
 */
public class SamClipIndelFraction extends AbstractCommandLineProgram
	{
	private enum Type {
		leftclip,
		rightclip,
		allclip,
		insert,
		deletion
		};
	private Type type=Type.allclip;
		
	private SamClipIndelFraction()
		{
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Extract clipping/indel fraction from BAM";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/SamClipIndelFraction";
		}

	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -t (type) one of: "+Arrays.toString(Type.values())+" (default:"+this.type+")");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
	
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"t:"))!=-1)
			{
			switch(c)
				{
				case 't': type= Type.valueOf(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		SamReader sfr=null;
		SAMRecordIterator iter=null;
		try
			{
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			if(opt.getOptInd()==args.length)
				{
				sfr=srf.open(SamInputResource.of(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{	
				String filename=args[opt.getOptInd()];
				sfr=srf.open(new File(filename));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			long total_bases_count=0L;
			long count_clipped_reads=0L;
			long count_clipped_left_reads=0L;
			long count_clipped_right_reads=0L;
			long count_unclipped_reads=0L;
			long count_unmapped_reads=0L;
			SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(sfr.getFileHeader().getSequenceDictionary());
			Counter<Integer> counter=new Counter<>();
			iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord record=iter.next();
				progress.watch(record);
				if(record.getReadUnmappedFlag())
					{
					++count_unmapped_reads;
					continue;
					}
				Cigar cigar=record.getCigar();
				int left_clip_length=0;
				int right_clip_length=0;
				int deletion_N_length=0;
				int deletion_D_length=0;
				int insertion_length=0;
				
				boolean onLeft=true;
				
				for(int i=0;i<  cigar.numCigarElements();++i)
					{
					CigarElement ce= cigar.getCigarElement(i);
					CigarOperator op = ce.getOperator();
					
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
					default: error("Bad type: "+type);return -1;
					}
				}
			progress.finish();
			
			System.out.println("##UNMAPPED_READS="+count_unmapped_reads);
			System.out.println("##MAPPED_READS="+(count_clipped_reads+count_unclipped_reads));
			System.out.println("##CLIPPED_READS="+count_clipped_reads);
			System.out.println("##CLIPPED_READS_5_PRIME="+count_clipped_left_reads);
			System.out.println("##CLIPPED_READS_3_PRIME="+count_clipped_right_reads);
			System.out.println("##UNCLIPPED_READS="+count_unclipped_reads);
			System.out.println("##COUNT_BASES="+total_bases_count);
			System.out.print("#");
			switch(type)
				{
				case leftclip:  System.out.print("CLIP_5_PRIME"); break;
				case rightclip:  System.out.print("CLIP_3_PRIME"); break;
				case allclip:  System.out.print("CLIP"); break;
				case deletion:  System.out.print("DELETION"); break;
				case insert:  System.out.print("INSERTION"); break;
				default: error("Bad type: "+type);return -1;
				}
			System.out.println("\tCOUNT\tFRACTION_OF_MAPPED_READS");
			
			for(Integer size: new TreeSet<Integer>(counter.keySet()))
				{
				System.out.print(size);
				System.out.print('\t');
				System.out.print(counter.count(size));
				System.out.print('\t');
				System.out.println(counter.count(size)/(double)(count_unclipped_reads+count_unclipped_reads));
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamClipIndelFraction().instanceMainWithExit(args);

	}

}
