package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import net.sf.picard.fastq.BasicFastqWriter;
import net.sf.picard.fastq.FastqConstants;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqGrep
	extends AbstractCommandLineProgram
	{
	private boolean inverse=false;
	private Map<String,Integer> readNames=new HashMap<String,Integer>(); 
	private int n_before_remove=-1;

	
	private FastqGrep()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqGrep";
		}
	
	@Override
	public String getProgramDescription() {
		return "Grep reads names in fastq";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (file) file containing a list of read names..");
		out.println(" -R (name) add the read.");
		out.println(" -n (int) when found, remove the read from the list of names when found more that 'n' time (increase speed)");
		out.println(" -V  invert");
		out.println(" -o (filename) output file. default: stdout.");
		super.printOptions(out);
		}
	
	private String getReadName(FastqRecord r)
		{
		return getReadName(r.getReadHeader());
		}
	
	private String getReadName(String s)
		{
		int beg=(s.startsWith(FastqConstants.SEQUENCE_HEADER)?1:0);
		int end=s.indexOf(' ');
		if(end==-1) end=s.length();
		s= s.substring(beg, end);
		return s;
		}
	private void run(FastqReader r,FastqWriter out)
		{
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			FastqRecord fastq=r.next();
			boolean keep=false;
			String readName=getReadName(fastq);
			Integer count=readNames.get(readName);
			if(count!=null)
				{
				keep=true;
				}
			if(inverse) keep=!keep;
			if(keep)
				{
				++nRec;
				out.write(fastq);
				}
			
			if(n_before_remove!=-1 && !inverse && keep)
				{
				count++;
				if(count>=n_before_remove)
					{
					readNames.remove(readName);
					if(readNames.isEmpty()) break;
					}
				else
					{
					readNames.put(readName,count);
					}
				}
			
			
			}
		info("Done. N-Reads:"+nRec);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File fileout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:f:R:n"))!=-1)
			{
			switch(c)
				{
				case 'n': n_before_remove=Integer.parseInt(opt.getOptArg()); break;
				case 'V': inverse=true;break;
				case 'R': readNames.put(getReadName(opt.getOptArg()),0);break;
				case 'f':
					{
					BufferedReader in=null;
					try
						{
						in=IOUtils.openURIForBufferedReading(opt.getOptArg());
				    	String line;
				    	while((line=in.readLine())!=null)
				    		{
				    		line=line.trim();
				    		if(line.isEmpty()) continue;
				    		readNames.put(getReadName(line),0);
				    		}
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					finally
						{
						CloserUtil.close(in);
						}
					break;
					}
				case 'o': fileout=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(readNames.isEmpty())
    		{
    		warning("no read name found.");
    		}
		
		FastqWriter out=null;
		try
			{
			if(fileout!=null)
				{
				info("Writing to "+fileout);
				out=new BasicFastqWriter(fileout);
				}
			else
				{
				info("Writing to stdout");
				out=new BasicFastqWriter(System.out);
				}
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				FastqReader fqR=new FourLinesFastqReader(System.in);
				run(fqR,out);
				fqR.close();
				}
			else for(int optind=opt.getOptInd(); optind < args.length; ++optind)
				{
				File f=new File(args[optind]);
				info("Reading from "+f);
				FastqReader fqR=new FourLinesFastqReader(f);
				run(fqR,out);
				fqR.close();
				}
			CloserUtil.close(out);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	public static void main(String[] args) {
		new FastqGrep().instanceMainWithExit(args);

	}

}
