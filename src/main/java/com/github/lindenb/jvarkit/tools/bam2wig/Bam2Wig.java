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
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class Bam2Wig extends AbstractCommandLineProgram
	{
	private int WINDOW_SIZE=100;
	private int WINDOW_SHIFT=25;
	private int min_qual=0;
	private boolean custom_track=false;
	private boolean cast_to_integer=false;
	private int min_gap=200;
	private int min_depth=0;
	private PrintWriter pw;
	
	private Bam2Wig()
		{
		
		}
	
	@Override
	public String getProgramDescription() {
		return "Bam to fixedStep Wiggle converter. Parses the cigar String to get the depth."+
				"Memory intensive: must alloc sizeof(int)*size(chrom)";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/Bam2Wig";
		}
	
	private void run(SamReader sfr)
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
				if(rec.getMappingQuality()< min_qual) continue;
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
				info("Allocating int["+prev_ssr.getSequenceLength()+"]");
				array=new int[prev_ssr.getSequenceLength()];
				info("Allocating : Done.");
				Arrays.fill(array, 0);
				}
			
			Cigar cigar=rec.getCigar();
			if(cigar==null) continue;
    		int refpos1=rec.getAlignmentStart();
    		for(CigarElement ce:cigar.getCigarElements())
    			{
				switch(ce.getOperator())
					{
					case H:break;
					case S:break;
					case I:break;
					case P:break;
					case N:// reference skip
					case D://deletion in reference
						{
    					refpos1+=ce.getLength();
						break;
						}
					case M:
					case EQ:
					case X:
						{
						for(int i=0;i< ce.getLength() && refpos1<= array.length;++i)
    		    			{
							if(refpos1>= 1 && refpos1<=array.length)
								{
								array[refpos1-1]++;
								}
    						refpos1++;
		    				}
						break;
						}
					default: throw new IllegalStateException(
							"Doesn't know how to handle cigar operator:"+ce.getOperator()+
							" cigar:"+cigar
							);

					}
    				
    			}
			

			
			
			}
		progess.finish();
		iter.close();
		pw.flush();
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -t print a UCSC custom track header.");
		out.println(" -s (int) window shift . Default:"+WINDOW_SHIFT);
		out.println(" -w (int) window size . Default:"+WINDOW_SIZE);
		out.println(" -q (int) min qual . Default:"+min_qual);
		out.println(" -i cast to integer.");
		out.println(" -g (int) minimal zero-coverage length before writing a new header. Default:"+min_gap);
		out.println(" -d (int) minimal depth before setting depth to zero . Default:"+min_depth);
		super.printOptions(out);	
		}
	
	@Override
	public int doWork(String[] args)
		{
	    GetOpt opt=new GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "ts:w:q:g:d:i"))!=-1)
			{
			switch(c)
				{
				case 'i': cast_to_integer=true;break;
				case 't':custom_track=true;break;
				case 's': WINDOW_SHIFT=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'w': WINDOW_SIZE=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'q': min_qual=Math.max(0, Integer.parseInt(opt.getOptArg())); break;
				case 'g': min_gap=Math.max(1, Integer.parseInt(opt.getOptArg())); break;
				case 'd': min_depth=Math.max(0, Integer.parseInt(opt.getOptArg())); break;
				default: 
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
				}
			}
		
		
		
	    
		SamReaderFactory sfrf= SamReaderFactory.makeDefault();
		sfrf.validationStringency( ValidationStringency.SILENT);
		SamReader in=null;
		this.pw=new PrintWriter(System.out);
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				in=sfrf.open(SamInputResource.of(System.in));
				run(in);
				in.close();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				in=sfrf.open(SamInputResource.of(new File(filename)));
				run(in);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			pw.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(pw);
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
