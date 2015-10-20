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
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Collection;

import htsjdk.variant.vcf.VCFCodec;
import com.github.lindenb.jvarkit.util.command.Command;
import htsjdk.samtools.util.CloserUtil;
import com.github.lindenb.jvarkit.io.IOUtils;

public class NoEmptyVCF extends AbstractNoEmptyVCF
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(NoEmptyVCF.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractNoEmptyVCF.AbstractNoEmptyVCFCommand
	 	{
    private void write(File f) throws IOException
    	{
    	LOG.info("fixing "+f);
		PrintWriter w=new PrintWriter(IOUtils.openFileForBufferedWriting(f));
		write(w);
		w.close();
    	}
    
    private void write(PrintWriter out) throws IOException
    	{
		out.println(VCFCodec.VCF4_MAGIC_HEADER);
		out.println("##source=noEmptyVCF");
		out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(!this.SAMPLES.isEmpty())
			{
			out.print("\tFORMAT");
			for(String S:this.SAMPLES)
				{
				out.print("\t");
				out.print(S);
				}
			}
		out.println();
		out.flush();
    	}
    @Override
    protected Collection<Throwable> call(String inputName) throws Exception {
    		PrintStream out= null;
			try
				{
				if(inputName!=null)
					{
					final File IN =new File(inputName);
					if(IN.exists())
						{
						LOG.info("opening "+IN+" ... ");
						Reader reader=IOUtils.openFileForBufferedReading(IN);
						int c=reader.read();
						reader.close();
						reader=null;
						if(c!=-1 && c!='#')
							{
							return wrapException("File "+IN+" doesn't start with # but with ascii("+c+")");
							}
						if(c==-1)
							{
							write(IN);
							}
						}
					else
						{
						LOG.error("File "+IN+" doesn't exists");
						write(IN);
						}
					}
				else
					{
					out = openFileOrStdoutAsPrintStream();
					LOG.info("reading from stdin ... ");
					int c= stdin().read();
					if(c!=-1 && c!='#')
						{
						LOG.warn("VCF doesn't start with # but with ascii("+c+")");
						out.print((char)c);
						IOUtils.copyTo(stdin(),out);
						return RETURN_OK;
						}
					if(c==-1)
						{
						LOG.warn("writing empty VCF");
						write(new PrintWriter(out));
						}
					else
						{
						out.print((char)c);
						IOUtils.copyTo(stdin(), out);
						}
					out.flush();
					}
				return RETURN_OK;
				}
			catch(IOException err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(out);
				}
			
			}
	 	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new NoEmptyVCF().instanceMainWithExit(args);
	}

}
