package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqToFasta
	extends AbstractCommandLineProgram
	{
	private int fastaLineLen=50;
	private boolean trim_after_space=false;
	
	private FastqToFasta()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqToFasta";
		}
	
	@Override
	public String getProgramDescription() {
		return "FastqToFasta";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (fileout) Filename output . Optional ");
		out.println(" -N (fasta line length)  Optional. Default: "+fastaLineLen);
		out.println(" -b trim fasta header after space.");
		super.printOptions(out);
		}
	
	
	
	private void run(FastqReader r,PrintStream out)
		{
		int wsp=0;
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			if(++nRec%1E6==0)
				{
				info("N-Reads:"+nRec);
				}
			FastqRecord fastq=r.next();
			out.print(">");
			if(!trim_after_space || (wsp=fastq.getReadHeader().indexOf(' '))==-1)
				{
				out.println(fastq.getReadHeader());
				}
			else
				{
				out.println(fastq.getReadHeader().substring(0, wsp));
				}
			
			int readLen=fastq.getReadString().length();
			int i=0;
			while(i< readLen)
				{
				int end=Math.min(i+fastaLineLen,readLen);
				out.println(fastq.getReadString().substring(i, end));
				i=end;
				}
			
			if(out.checkError()) break;
			}
		out.flush();
		info("Done. N-Reads:"+nRec);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File fileout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args, getGetOptDefault()+"o:N:b"))!=-1)
			{
			switch(c)
				{
				case 'o': fileout=new File(getopt.getOptArg());break;
				case 'b': trim_after_space=true;break;
				case 'N': fastaLineLen=Math.max(1,Integer.parseInt(getopt.getOptArg()));break;
				default: 
					{
					switch(handleOtherOptions(c, getopt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		PrintStream out=System.out;
		try
			{
			if(fileout!=null)
				{
				info("Writing to "+fileout);
				out=new PrintStream(IOUtils.openFileForWriting(fileout));
				}
			else
				{
				info("Writing to stdout");
				out=System.out;
				}
			
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				FastqReader fqR=new FourLinesFastqReader(System.in);
				run(fqR,out);
				fqR.close();
				}
			else for(int optind=getopt.getOptInd(); optind < args.length; ++optind)
				{
				File f=new File(args[optind]);
				info("Reading from "+f);
				FastqReader fqR=new FourLinesFastqReader(f);
				run(fqR,out);
				fqR.close();
				}
			out.flush();
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
		return 0;
		}
	
	public static void main(String[] args) {
		new FastqToFasta().instanceMainWithExit(args);

	}

}
