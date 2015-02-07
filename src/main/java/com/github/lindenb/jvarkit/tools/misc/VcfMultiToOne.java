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
import java.util.Collections;
import java.util.List;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/*
 * VcfMultiToOne
 */
public class VcfMultiToOne extends AbstractVCFFilter2
	{
	private boolean keep_no_call=true;
	private boolean keep_hom_ref=true;
	private boolean keep_non_available=true;
	
	private VcfMultiToOne()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfMultiToOne";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final String VCF_SAMPLE_NAME="SAMPLE";
		final String INFO_TAG="SAMPLENAME";
		VCFHeader srcHeader=in.getHeader();
		List<String> sampleNames=srcHeader.getSampleNamesInOrder();
		VCFHeader header=new VCFHeader(srcHeader.getMetaDataInInputOrder(),
				Collections.singleton(VCF_SAMPLE_NAME));
		header.addMetaDataLine(new VCFInfoHeaderLine(
				INFO_TAG,1,VCFHeaderLineType.String,"Sample Name from multi-sample vcf"
				));

		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		
		for(String sample:sampleNames)
			{
			header.addMetaDataLine(
				new VCFHeaderLine(getClass().getSimpleName()+".Sample",
				sample));
			}
		
		out.writeHeader(header);
		while(in.hasNext() && !System.out.checkError())
			{
			VariantContext ctx = in.next();
			if(sampleNames.isEmpty())
				{
				out.add(ctx);
				continue;
				}
			for(String sample:sampleNames)
				{
				Genotype g= ctx.getGenotype(sample);
				if(!g.isCalled() && !keep_no_call) continue;
				if(!g.isAvailable() && !keep_non_available) continue;
				if(g.isHomRef() && !keep_hom_ref) continue;
				
				
				GenotypeBuilder gb=new GenotypeBuilder(g);
				gb.name(VCF_SAMPLE_NAME);
				
				
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(INFO_TAG, sample);
				vcb.genotypes(gb.make());
				out.add(vcb.make());
				}
			}
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -c discard if variant is no-call");
		out.println(" -r discard if variant is hom-ref");
		out.println(" -a discard if variant is non-available");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "cra"))!=-1)
			{
			switch(c)
				{
				case 'c': this.keep_no_call=false; break;
				case 'r': this.keep_hom_ref=false; break;
				case 'a': this.keep_non_available=false; break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
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
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
