package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintStream;

import net.sf.picard.fastq.FastqConstants;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;

public class FastqRevComp extends AbstractCommandLineProgram
	{

	private FastqRevComp()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqRevComp";
		}
	
	@Override
	public String getProgramDescription() {
		return "produces a reverse-complement fastq (for mate pair alignment see http://seqanswers.com/forums/showthread.php?t=5085 )";
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
		while(r.hasNext())
			{
			if(++nRec%1E6==0)
				{
				info("N-Reads:"+nRec);
				}
			FastqRecord fastq=r.next();
			
			out.print(FastqConstants.SEQUENCE_HEADER);
			out.println(fastq.getReadHeader());
			String s=fastq.getReadString();
			for(int i=s.length()-1;i>=0;i--)
				{
				out.print(AcidNucleics.complement(s.charAt(i)));
				}
			out.println();
			
			out.print(FastqConstants.QUALITY_HEADER);
			s=fastq.getBaseQualityHeader();
			if(s!=null) out.print(s);
			out.println();
			
			s=fastq.getBaseQualityString();
			
			for(int i=s.length()-1;i>=0;i--)
				{
				out.print(s.charAt(i));
				}
			out.println();
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
				FastqReader fqR=new FastqReader(new BufferedReader(new InputStreamReader(System.in)));
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
		new FastqRevComp().instanceMainWithExit(args);

	}

}
