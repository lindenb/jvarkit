package com.github.lindenb.jvarkit.tools.fastq;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;

import net.sf.picard.fastq.FastqConstants;
import net.sf.picard.fastq.FastqRecord;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class ConvertPhred64toFastq33 extends AbstractCommandLineProgram
	{
	private PrintStream pw= System.out;
	private ConvertPhred64toFastq33()
		{
		
		}
	@Override
	public String getProgramDescription() {
		return "Convert Illumina Fastq 64 encoding to Fastq 33.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	private void convert(InputStream in) throws IOException
		{
		FastqReader r=new FourLinesFastqReader(in);
		while(r.hasNext() && !pw.checkError())
			{
			FastqRecord rec=r.next();
			byte quals[]=rec.getBaseQualityString().getBytes();
			for(int i=0;i< quals.length;++i )
				{
				quals[i]=(byte)(quals[i]-64+33);
				if(quals[i]<33 || quals[i]>126)
					{
					r.close();
					throw new IOException("q="+(int)quals[i]);
					}
				}
			String name=rec.getReadHeader();
			int diez=name.indexOf('#');
			if(diez!=-1) name=name.substring(0, diez);
	        pw.print(FastqConstants.SEQUENCE_HEADER);
	        pw.println(name);
	        pw.println(rec.getReadString());
	        pw.print(FastqConstants.QUALITY_HEADER);
	        pw.println(rec.getBaseQualityHeader() == null || rec.getReadHeader().equals(rec.getBaseQualityHeader())? "" : rec.getBaseQualityHeader());
	        pw.println(new String(quals));
			}
		r.close();
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
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
		
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				convert(System.in);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					InputStream in=IOUtils.openURIForReading(filename);
					convert(in);
					in.close();
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
			
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new ConvertPhred64toFastq33().instanceMainWithExit(args);
	}

}
