package com.github.lindenb.jvarkit.tools.fastq;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.fastq.FastqRecord;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqEntropy extends AbstractCommandLineProgram
	{
	private PrintStream pw= System.out;
	private Counter<Long> length2count=new Counter<Long>();
	
	private static class BestCompressionOutputStream extends GZIPOutputStream
		{
		BestCompressionOutputStream() throws IOException
			{
			super(new NullOuputStream());
			def.setLevel(Deflater.BEST_COMPRESSION);
			
			}
		public long getByteWrittenCount()
			{
			return NullOuputStream.class.cast(super.out).getByteWrittenCount();
			}
		}
	private FastqEntropy()
		{
		
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqEntropy";
		}
	@Override
	public String getProgramDescription() {
		return "Compute the Entropy of a Fastq file (distribution of the length(gzipped(sequence)))";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	private void convert(InputStream in) throws IOException
		{
		FastqReader r=new FourLinesFastqReader(in);
		while(r.hasNext())
			{
			FastqRecord rec=r.next();
			BestCompressionOutputStream gzout=new BestCompressionOutputStream();
			gzout.write(rec.getBaseQualityString().getBytes());
			gzout.flush();
			gzout.close();
			this.length2count.incr(gzout.getByteWrittenCount());
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
			for(Long n:this.length2count.keySetIncreasing())
				{
				pw.print(n);
				pw.print('\t');
				pw.println(this.length2count.count(n));
				if(pw.checkError()) break;
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
		new FastqEntropy().instanceMainWithExit(args);
	}

}
