/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum PhD.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.List;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqRevComp extends AbstractFastqRevComp
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FastqRevComp.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractFastqRevComp.AbstractFastqRevCompCommand
	 	{
		
		private void run(FastqReader r,PrintStream out)
			{
			String s;
			long nRec=0L;
			r.setValidationStringency(ValidationStringency.LENIENT);
			while(r.hasNext())
				{
				if(++nRec%1E6==0)
					{
					LOG.info("N-Reads:"+nRec);
					}
				FastqRecord fastq=r.next();
				
				
				out.print(FastqConstants.SEQUENCE_HEADER);
				out.println(fastq.getReadHeader());
				s=fastq.getReadString();
	
				if((this.only_R2 && nRec%2==1) || (this.only_R1 && nRec%2==0) ) //interleaced
					{
					out.print(s);
					}
				else
					{
					for(int i=s.length()-1;i>=0;i--)
						{
						out.print(AcidNucleics.complement(s.charAt(i)));
						}
					}
				out.println();
				
				out.print(FastqConstants.QUALITY_HEADER);
				s=fastq.getBaseQualityHeader();
				if(s!=null) out.print(s);
				out.println();
				
				s=fastq.getBaseQualityString();
				if((this.only_R2 && nRec%2==1) || (this.only_R1 && nRec%2==0) ) //interleaced
					{
					out.print(s);
					}
				else
					{
					for(int i=s.length()-1;i>=0;i--)
						{
						out.print(s.charAt(i));
						}
					}
				out.println();
				if(out.checkError()) break;
				}
			out.flush();
			LOG.info("Done. N-Reads:"+nRec);
			}
		
		@Override
		public Collection<Throwable> call() throws Exception
			{
			if(only_R1 && only_R2)
				{
				return wrapException("Both options -1 && -2 used.");
				}
			final List<String> args = getInputFiles();
			PrintStream out=System.out;
			try
				{
				out = super.openFileOrStdoutAsPrintStream();
				
				if(args.isEmpty())
					{
					LOG.info("Reading from stdin");
					FastqReader fqR=new FourLinesFastqReader(System.in);
					run(fqR,out);
					fqR.close();
					}
				else for(String arg: args)
					{
					File f=new File(arg);
					LOG.info("Reading from "+f);
					FastqReader fqR=new FourLinesFastqReader(f);
					run(fqR,out);
					fqR.close();
					}
				out.flush();
				return RETURN_OK;
				}
			catch(Exception err)
				{
				LOG.error(err);
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(out);
				}
			}
	 	}
	
	public static void main(String[] args)
		{
		new FastqRevComp().instanceMainWithExit(args);
		}
	}
