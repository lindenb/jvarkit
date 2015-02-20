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
* 2015 : moving to knime
* 2014 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfTail
	extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	/* number of row to limit */
	private int count=10;
	/** output file : stdout if null */
	private File outputFile=null;
	/** number of variants filterer */
	private int countFilteredVariants=0;

	
	public VcfTail()
		{
		}
	
	public void setCount(int count) {
		this.count = Math.max(0,count);
		}	

	@Override
	public String getProgramDescription() {
		return "Print last lines of a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfTail";
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		VcfIterator vcfIn=null;
		try
			{
			if(args.isEmpty())
				{
				vcfIn = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				vcfIn= VCFUtils.createVcfIterator(args.get(0));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			this.filterVcfIterator(vcfIn);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			}
		
		}

	public int geVariantCount()
		{
		return this.countFilteredVariants;
		}
	
	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}

	public File getOutputFile() {
		return outputFile;
		}

	private void filterVcfIterator(VcfIterator in) throws IOException
		{
		this.countFilteredVariants=0;
		VariantContextWriter out = null;
		try {
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			
			if(getOutputFile()==null)
				{
				out = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+getOutputFile());
				out = VCFUtils.createVariantContextWriter(getOutputFile());
				}
			
			out.writeHeader(h2);
			LinkedList<VariantContext> L=new LinkedList<VariantContext>();
			while(in.hasNext() && L.size()< this.count)
				{	
				L.add(progess.watch(in.next()));
				}
			while(in.hasNext())
				{
				L.add(progess.watch(in.next()));
				L.removeFirst();
				}
			for(VariantContext ctx:L)
				{
				out.add(ctx);
				}
			progess.finish();
			this.countFilteredVariants=L.size();
			}
		finally
			{
			CloserUtil.close(out);
			out=null;
			}
		}
		
		@Override
		public void printOptions(PrintStream out)
			{
			out.println(" -n (int) output size. Optional. Default:"+this.count);
			out.println(" -o (fileout). Optional. Default: stdout");
			super.printOptions(out);
			}
		
		@Override
		public int doWork(String[] args)
			{
			com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
			int c;
			while((c=opt.getopt(args,getGetOptDefault()+ "n:o:"))!=-1)
				{
				switch(c)
					{
					case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
					case 'n': this.setCount(Integer.parseInt(opt.getOptArg())); break;
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
			
			this.initializeKnime();
			List<String> L=new ArrayList<String>();
			for(int i=opt.getOptInd();i<args.length;++i)
				{
				L.add(args[i]);
				}
			return this.executeKnime(L);
			}
		@Override
		public int initializeKnime() {
			return 0;
			}
		@Override
		public void disposeKnime() {
			}
		
		@Override
		public void checkKnimeCancelled() {
			
			}
		
		
		
		public static void main(String[] args)
			{
			new VcfTail().instanceMainWithExit(args);
			}
		}
