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
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class PadEmptyFastq extends AbstractCommandLineProgram
	{
	private static final int DEFAULT_LENGTH=50;
	private int N=-1;//default , will use the first read length
	
	
	private PadEmptyFastq()
		{
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
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/PadEmptyFastq";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Pad empty fastq sequence/qual with N/#";
		}

	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (Filename out). Default: stdout");
		out.println(" -N (integer). number of bases/qual to be added.  Default: length of the first read  or "+DEFAULT_LENGTH);
		super.printOptions(out);
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		
		File fileOut=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:N:"))!=-1)
			{
			switch(c)
				{
				case 'o':fileOut=new File(opt.getOptArg());break;
				case 'N':N=Math.max(0,Integer.parseInt(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					break;
					}
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
				FastqReader fqr=new FourLinesFastqReader(System.in);
				copyTo(fqr,fqw);
				fqr.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					FastqReader fqr=new FourLinesFastqReader(new File(filename));
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
