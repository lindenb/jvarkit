/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;


@Program(name="fastqphred64to33",
	keywords={"fastq"},
	description="Convert Illumina Fastq 64 encoding to Fastq 33"
	)
public class ConvertPhred64toFastq33 extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertPhred64toFastq33.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	
	private PrintStream pw= System.out;
	private ConvertPhred64toFastq33()
		{
		
		}

	
	private void convert(InputStream in) throws IOException
		{
		FastqReader r=new FourLinesFastqReader(in);
		while(r.hasNext() && !pw.checkError())
			{
			final FastqRecord rec=r.next();
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
			String name=rec.getReadName();
			int diez=name.indexOf('#');
			if(diez!=-1) name=name.substring(0, diez);
	        pw.print(FastqConstants.SEQUENCE_HEADER);
	        pw.println(name);
	        pw.println(rec.getReadString());
	        pw.print(FastqConstants.QUALITY_HEADER);
	        pw.println(rec.getBaseQualityHeader() == null || rec.getReadName().equals(rec.getBaseQualityHeader())? "" : rec.getBaseQualityHeader());
	        pw.println(new String(quals));
			}
		r.close();
		}
	@Override
	public int doWork(List<String> args) {		
		try
			{
			this.pw = super.openFileOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				convert(stdin());
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					InputStream in=IOUtils.openURIForReading(filename);
					convert(in);
					in.close();
					}
				}
			this.pw.flush();
			this.pw.close();
			this.pw = null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new ConvertPhred64toFastq33().instanceMainWithExit(args);
	}

}
