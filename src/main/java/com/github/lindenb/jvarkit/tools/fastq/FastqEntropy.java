package com.github.lindenb.jvarkit.tools.fastq;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.fastq.FastqRecord;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.semontology.Term;

@Program(name="fastqentropy",
	description="Compute the Entropy of a Fastq file (distribution of the length(gzipped(sequence))",
	keywords={"fastq"},
	terms={Term.ID_0000005}
	)
public class FastqEntropy extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqEntropy.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File fileout = null;

	private PrintStream pw= System.out;
	private final Counter<Long> length2count=new Counter<Long>();
	
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
	public int doWork(List<String> args) {
		try
			{
			this.pw = super.openFileOrStdoutAsPrintStream(this.fileout);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				convert(stdin());
				}
			else
				{
				for(String filename: args)
					{
					LOG.info("Reading from "+filename);
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
			pw.flush();
			pw.close();
			pw=null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
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
