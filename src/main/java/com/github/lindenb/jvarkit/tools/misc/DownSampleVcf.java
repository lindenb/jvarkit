/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class DownSampleVcf extends AbstractVCFFilter2
	{
	private int reservoir_size=10;
	private long seed=System.currentTimeMillis();
	private DownSampleVcf()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "DownSample a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/DownSampleVcf";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		Random rand=new Random(this.seed);
		List<VariantContext>  buffer=new ArrayList<VariantContext>(this.reservoir_size);
		VCFHeader h2=new VCFHeader(in.getHeader());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(in.getHeader());
		out.writeHeader(h2);
		if(this.reservoir_size!=0)
			{
			while(in.hasNext())
				{	
				if(buffer.size() < this.reservoir_size)
					{
					buffer.add(progess.watch(in.next()));
					}
				else
					{
					buffer.set(rand.nextInt(buffer.size()), progess.watch(in.next()));
					}
				}
			
			}
		for(VariantContext ctx:buffer)
			{
			out.add(ctx);
			}
		progess.finish();
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -N (long) random seed. Optional.");
		out.println(" -n (int) output size. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "N:n:"))!=-1)
			{
			switch(c)
				{
				case 'N': seed=Long.parseLong(opt.getOptArg()); break;
				case 'n': reservoir_size=Math.max(0, Integer.parseInt(opt.getOptArg())); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		return doWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new DownSampleVcf().instanceMainWithExit(args);
		}
	}
