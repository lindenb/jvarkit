package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;

import net.sf.picard.fastq.BasicFastqWriter;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;

public class PadEmptyFastq extends AbstractCommandLineProgram
	{
	private static final int DEFAULT_LENGTH=50;
	private int N=-1;//default , will use the first read length
	
	
	private PadEmptyFastq()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/PadEmptyFastq";
	}
	
	@Override
	public String getProgramDescription()
		{
		return "Pad empty fastq sequence/qual with N/#";
		}
	
	private void copyTo(FastqReader r,FastqWriter w)
		{
		int padLength=this.N;
		long nReads=0L;
		long nFill=0L;
		String fillN=null;
		String fillQ=null;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			FastqRecord rec=r.next();
			
			
			if(++nReads%1E6==0)
				{
				info("Read "+nReads +" reads. empty reads="+nFill);
				}
			if(rec.getReadString().isEmpty())
				{
				++nFill;
				if(padLength<1)
					{
					padLength=DEFAULT_LENGTH;
					}
				if(fillN==null)
					{
					StringBuilder b1=new StringBuilder();
					while(b1.length()< padLength) b1.append("N");
					fillN=b1.toString();
					fillQ=fillN.replace('N', '#');
					}
				
				rec=new FastqRecord(
						rec.getReadHeader(),
						fillN,
						rec.getBaseQualityHeader(),
						fillQ
						);
				}
			else if(padLength<1)
				{
				padLength=rec.getReadString().length();
				}
			w.write(rec);
			}
		info("Done. Read "+nReads +" reads. empty reads="+nFill);
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -h get help (this screen)");
		out.println(" -o (Filename out). Default: stdout");
		out.println(" -N (integer). number of bases/qual to be added.  Default: length of the first read  or "+DEFAULT_LENGTH);

		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		File fileOut=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, "hvL:o:N:"))!=-1)
			{
			switch(c)
				{
				case 'o':fileOut=new File(opt.getOptArg());break;
				case 'N':N=Math.max(0,Integer.parseInt(opt.getOptArg()));break;
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(opt.getOptArg()));break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+opt.getOptOpt());return -1;
				}
			}
		
		FastqWriter fqw=null;		
		try
			{
			
			if(fileOut==null)
				{
				info("writing to stdout");
				fqw=new BasicFastqWriter(System.out);
				}
			else
				{
				info("writing to "+fileOut);
				fqw=new FastqWriterFactory().newWriter(fileOut);
				}
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				FastqReader fqr=new FastqReader(System.in);
				copyTo(fqr,fqw);
				fqr.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					FastqReader fqr=new FastqReader(new File(filename));
					copyTo(fqr,fqw);
					fqr.close();
					}
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
			CloserUtil.close(fqw);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new PadEmptyFastq().instanceMainWithExit(args);
		}

}
