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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfMultiToOneAllele
	extends AbstractVCFFilter3
	{

	
	 public VcfMultiToOneAllele()
		{
		}
	 
	@Override
	public String getProgramDescription() {
		return " 'one variant with N ALT alleles' to 'N variants with one ALT'";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfMultiToOneAllele";
		}
	
	
	@Override
	protected void doWork(
				String inpuSource,
				VcfIterator in,
				VariantContextWriter out
				)
			throws IOException {
		final String TAG="VCF_MULTIALLELIC_SRC";
		final List<String> noSamples=Collections.emptyList();
		final List<Genotype> noGenotypes=Collections.emptyList();
		VCFHeader header=in.getHeader();
		Set<VCFHeaderLine> metaData=new HashSet<>(header.getMetaDataInInputOrder());
		
		
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		metaData.add(new VCFInfoHeaderLine(TAG, 1, VCFHeaderLineType.String, "The variant was processed with VcfMultiAlleleToOneAllele and contained the following alleles."));
		
		VCFHeader h2=new VCFHeader(
				metaData,
				noSamples
				);
		
		SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		out.writeHeader(h2);
		while(in.hasNext())
			{
			VariantContext ctx=progess.watch(in.next());
			List<Allele> alleles = ctx.getAlternateAlleles();
			if(alleles.isEmpty())
				{
				warning("Remove no ALT variant:"+ctx);
				continue;
				}
			else if(alleles.size()==1)
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.genotypes(noGenotypes);
				out.add(vcb.make());
				incrVariantCount();
				}
			else
				{
				final Map<String,Object> attributes = ctx.getAttributes();
				StringBuilder sb=new StringBuilder();
				for(int i=0;i< alleles.size();++i)
					{
					if(sb.length()>0) sb.append("|");
					sb.append(alleles.get(i).getDisplayString());
					}
				String altAsString= sb.toString();
				
				for(int i=0;i< alleles.size();++i)
					{
					Allele a=alleles.get(i);
					VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.alleles(Arrays.asList(ctx.getReference(),a));
					for(String attid:attributes.keySet())
						{
						VCFInfoHeaderLine info = header.getInfoHeaderLine(attid);
						if(info==null) throw new IOException("Cannot get header INFO tag="+attid);
						if(info.getCountType()!=VCFHeaderLineCount.A) continue;
						Object o = 	attributes.get(attid);
						if(!(o instanceof List)) throw new IOException("For INFO tag="+attid+" got "+o.getClass()+" instead of List");
						@SuppressWarnings("rawtypes")
						List list = (List)o;
						if(i>=list.size()) throw new IOException("For INFO tag="+alleles.size()+" got "+alleles.size()+" ALT, incompatible with "+list.toString());
						vcb.attribute(attid, list.get(i));
						vcb.attribute(TAG,altAsString);
						incrVariantCount();
						}
					out.add(vcb.make());
					}
				}
			}
		progess.finish();
		}

		
	@Override
	public void printOptions(PrintStream out)
		{
		super.printOptions(out);
		}
		
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
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
		
		return mainWork(opt.getOptInd(), args);
		}
	
	
	public static void main(String[] args)
		{
		new VcfMultiToOneAllele().instanceMainWithExit(args);
		}
	
	}
