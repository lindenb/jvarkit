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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

import htsjdk.variant.vcf.VCFCodec;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
 BEGIN_DOC
 
## History
 
 2017:  moved to jcommander
 
 END_DOC
 */
@Program(name="noemptyvcf",
	description="If VCF is empty or doesn't exists, create a dummy one",
	deprecatedMsg="Was developped at the time where VEP didn't send an output if there was no variant, just a header in the source vcf."
	)
public class NoEmptyVCF extends Launcher
	{
	private static final Logger LOG=Logger.build(NoEmptyVCF.class).make();
    
	@Parameter(names={"-s","--sample"},description="Add this genotyped samples")
    private List<String> SAMPLES=new ArrayList<String>();
	@Parameter(names={"-o","--out"},description="Output VCF or stdout")
    private File outFile=null;

	
    
    private void writeEmptyVcf() throws IOException
    	{
		final PrintWriter w;
		
		if(outFile!=null)
			{
			w = new PrintWriter(IOUtils.openFileForBufferedWriting(outFile));
			}
		else
			{
			w= new PrintWriter(stdout());
			}
    	
		w.println(VCFCodec.VCF4_MAGIC_HEADER);
		w.println("##source=noEmptyVCF");
		w.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(!this.SAMPLES.isEmpty())
			{
			w.print("\tFORMAT");
			for(final String S:this.SAMPLES)
				{
				w.print("\t");
				w.print(S);
				}
			}
		w.println();
		w.flush();
    	}
    @Override
    public int doWork(final List<String> args) {
    	String inputName=super.oneFileOrNull(args);;
		try
			{
			if(inputName!=null)
				{
				final File IN=new File(inputName);
				if(IN.exists())
					{
					LOG.info("opening "+IN+" ... ");
					Reader reader=IOUtils.openFileForBufferedReading(IN);
					int c=reader.read();
					reader.close();
					reader=null;
					if(c!=-1 && c!='#')
						{
						LOG.fatal("File "+IN+" doesn't start with # but with ascii("+c+")");
						return -1;
						}
					if(c==-1)
						{
						writeEmptyVcf();
						}
					}
				else
					{
					LOG.warn("File "+IN+" doesn't exists");
					writeEmptyVcf();
					}
				}
			else
				{
				LOG.info("reading from stdin ... ");
				int c=stdin().read();
				if(c!=-1 && c!='#')
					{
					LOG.warn("VCF doesn't start with # but with ascii("+c+")");
					System.out.print((char)c);
					IOUtils.copyTo(stdin(), stdout());
					return 0;
					}
				if(c==-1)
					{
					LOG.warn("writing empty VCF");
					writeEmptyVcf();
					}
				else
					{
					stdout().print((char)c);
					IOUtils.copyTo(stdin(), stdout());
					}
				}
			}
		catch(final IOException err)
			{
			LOG.fatal(err);
			return -1;
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(final String[] args) {
		new NoEmptyVCF().instanceMainWithExit(args);

	}

}
