/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.List;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqToFasta
	extends AbstractFastqToFasta
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FastqToFasta.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFastqToFasta.AbstractFastqToFastaCommand
	 	{		
	
	
	private void run(FastqReader r,PrintStream out)
		{
		int wsp=0;
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			if(++nRec%1E6==0)
				{
				LOG.info("N-Reads:"+nRec);
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
		LOG.info("Done. N-Reads:"+nRec);
		}
	
	@Override
	public Collection<Throwable> call() throws Exception
			{
			final List<String> args = getInputFiles();
			PrintStream out=null;
			FastqReader fqR = null;
			try
				{
				out =  openFileOrStdoutAsPrintStream();
				
				if(args.isEmpty())
					{
					LOG.info("Reading from stdin");
					 fqR=new FourLinesFastqReader(stdin());
					run(fqR,out);
					fqR.close();
					fqR =null;
					}
				else for(String arg: args)
					{
					File f=new File(arg);
					LOG.info("Reading from "+f);
					fqR=new FourLinesFastqReader(f);
					run(fqR,out);
					fqR.close();
					fqR =null;
					}
				out.flush();
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(out);
				CloserUtil.close(fqR);
				}
			}
	 	}
	
	public static void main(String[] args) {
		new FastqToFasta().instanceMainWithExit(args);

	}

}
