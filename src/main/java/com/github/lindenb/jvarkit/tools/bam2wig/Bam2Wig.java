package com.github.lindenb.jvarkit.tools.bam2wig;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;


import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;

public class Bam2Wig extends AbstractCommandLineProgram
	{
	private int WINDOW_SIZE=100;
	private int WINDOW_SHIFT=25;
	private int min_qual=0;
	private boolean custom_track=false;
	private boolean cast_to_integer=false;
	private int min_gap=200;
	private int min_depth=0;
	
	
	@Override
	public String getProgramDescription() {
		return "Bam to fixedStep Wiggle converter. Parses the cigar String to get the depth."+
				"Memory intensive: must alloc sizeof(int)*size(chrom)";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/Bam2Wig";
		}
	
	private void run(SAMFileReader sfr)
		{
		SAMSequenceRecord prev_ssr=null;
		sfr.setValidationStringency(ValidationStringency.LENIENT);
		SAMSequenceDictionary dict=sfr.getFileHeader().getSequenceDictionary();
		PrintWriter w=new  PrintWriter(System.out);
		SAMRecordIterator iter=sfr.iterator();
		int array[]=null;
		long nReads=0;
		if(custom_track)
			{
			w.println("track type=wiggle_0 name=\"__REPLACE_WIG_NAME__\" description=\"__REPLACE_WIG_DESC__\"");
			}
		
		for(;;)
			{
			SAMRecord rec=null;
			if(iter.hasNext())
				{
				rec=iter.next();
				if(++nReads%1E6==0)
					{
					info("nReads: "+nReads);
					}
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
		 								w.println(0);
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
			 					w.println(
				 						"fixedStep chrom="+prev_ssr.getSequenceName()+
				 						" start="+(start0+1)+
				 						" step="+WINDOW_SHIFT +" span="+ WINDOW_SIZE);
			 					}
		 					
		 					if(cast_to_integer)
		 						{
		 						w.println((int)(sum/n));
		 						}
		 					else
		 						{
		 						w.println((float)(sum/n));
		 						}
		 					
		 				if(w.checkError()) break;
		 				start0+=WINDOW_SHIFT;
						}
					array=null;
					System.gc();
					prev_ssr=null;
					}
				if(rec==null) break;
				if(w.checkError()) break;
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
		iter.close();
		w.flush();
		w.close();
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
	    GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args,getGetOptDefault()+ "ts:w:q:g:d:i"))!=-1)
			{
			switch(c)
				{
				case 'i': cast_to_integer=true;break;
				case 't':custom_track=true;break;
				case 's': WINDOW_SHIFT=Math.max(1, Integer.parseInt(getopt.getOptArg())); break;
				case 'w': WINDOW_SIZE=Math.max(1, Integer.parseInt(getopt.getOptArg())); break;
				case 'q': min_qual=Math.max(0, Integer.parseInt(getopt.getOptArg())); break;
				case 'g': min_gap=Math.max(1, Integer.parseInt(getopt.getOptArg())); break;
				case 'd': min_depth=Math.max(0, Integer.parseInt(getopt.getOptArg())); break;
				default: 
					switch(handleOtherOptions(c, getopt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
				}
			}
		
		
		
	    
		SAMFileReader samFileReader=null;
		try
			{
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				samFileReader=new SAMFileReader(System.in);
				}
			else if(getopt.getOptInd()+1==args.length)
				{
				File bamFile=new File(args[getopt.getOptInd()]);
				info("Reading from "+bamFile);
				samFileReader=new SAMFileReader(bamFile);
				}
			else
				{
				System.err.println("illegal number of arguments.");
				return -1;
				}
			
			run(samFileReader);
			
			
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samFileReader);
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Bam2Wig().instanceMainWithExit(args);
		}

	}
