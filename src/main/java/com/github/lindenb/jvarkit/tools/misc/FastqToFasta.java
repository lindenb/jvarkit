package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;

import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;

public class FastqToFasta extends AbstractCommandLineProgram
	{

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
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -o (fileout) Filename output . Optional ");
		}
	
	
	
	private void run(FastqReader r,PrintStream out)
		{
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
			out.println(fastq.getReadHeader());
			int readLen=fastq.getReadString().length();
			int i=0;
			while(i< readLen)
				{
				int end=Math.min(i+50,readLen);
				out.println(fastq.getReadString().substring(i, end));
				i=end;
				}
			if(i%50!=0) out.println();
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
		while((c=getopt.getopt(args, "hvL:o:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(getopt.getOptArg()));break;
				case 'o': fileout=new File(getopt.getOptArg());break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
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
				FastqReader fqR=new FastqReader(System.in);
				run(fqR,out);
				fqR.close();
				}
			else for(int optind=getopt.getOptInd(); optind < args.length; ++optind)
				{
				File f=new File(args[optind]);
				info("Reading from "+f);
				FastqReader fqR=new FastqReader(f);
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
