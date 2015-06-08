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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class Biostar130456 extends AbstractCommandLineProgram
	{
	private final static String SAMPLE_TAG="__SAMPLE__";

	private Biostar130456()
		{
		}
	@Override
	public String getProgramDescription()
		{
		return " Individual VCF files from main VCF file. See   https://www.biostars.org/p/130456";
		}
	
	@Override
	public String getProgramName()
		{
		return "Biostar130456";
		}
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar130456";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -p (pattern) output file pattern. Must contain the word "+SAMPLE_TAG);
		out.println(" -x remove uncalled genotypes");
		out.println(" -z remove homzygote REF/REF");
		super.printOptions(out);
		}
	

	
	
	
	@Override
	public int doWork(String[] args)
		{
		boolean remove_uncalled =false;
		boolean remove_homref =false;
		String filepattern=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"p:xz"))!=-1)
			{
			switch(c)
				{
				case 'p': filepattern=opt.getOptArg();break;
				case 'x': remove_uncalled=true;break;
				case 'z': remove_homref=true;break;
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
		if(filepattern==null || !filepattern.contains(SAMPLE_TAG))
			{
			error("File pattern is missing "+SAMPLE_TAG);
			return -1;
			}
		
		VcfIterator in=null;
		try
			{
			if(args.length==opt.getOptInd())
				{
				info("Reading from <stdin>");
				in = VCFUtils.createVcfIteratorStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				in = VCFUtils.createVcfIterator(filename);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			VCFHeader header=in.getHeader();
			Set<String> samples = new HashSet<String>(header.getSampleNamesInOrder());
			Map<String,VariantContextWriter> sample2writer=new HashMap<String,VariantContextWriter>(samples.size());

			if(samples.isEmpty())
				{
				info("VCF doesn't contain any sample");
				return -1;
				}
			info("N sample:"+samples.size());
			for(String sample:samples)
				{
				
				VCFHeader h2=new VCFHeader(
						header.getMetaDataInInputOrder(),
						Collections.singleton(sample)
						);
				h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
				h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
				h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
				h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
				String sampleFile= filepattern.replaceAll(SAMPLE_TAG,sample);
				System.out.println(sampleFile);
				File fout = new File(sampleFile);
				if(fout.getParentFile()!=null) fout.getParentFile().mkdirs();
				VariantContextWriter w= VCFUtils.createVariantContextWriter(fout);
				w.writeHeader(h2);
				sample2writer.put(sample, w);
				}
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			while(in.hasNext())
				{
				VariantContext ctx= progress.watch(in.next());
				for(String sample: samples)
					{
					Genotype g= ctx.getGenotype(sample);
					if(g==null) continue;
					if(remove_uncalled && (!g.isAvailable() || !g.isCalled() || g.isNoCall()))
						{
						continue;
						}
					if(remove_homref && g.isHomRef()) continue;
					VariantContextWriter w= sample2writer.get(sample);
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					GenotypeBuilder gb=new GenotypeBuilder(g);
					vcb.genotypes(Collections.singletonList(gb.make()));
					VariantContext ctx2= vcb.make();
					w.add(ctx2);
					}
				}
			for(String sample:samples)
				{
				info("Closing for sample "+sample);
				VariantContextWriter w= sample2writer.get(sample);
				w.close();
				}
			progress.finish();
			return 0;
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar130456().instanceMainWithExit(args);
		}

	}
