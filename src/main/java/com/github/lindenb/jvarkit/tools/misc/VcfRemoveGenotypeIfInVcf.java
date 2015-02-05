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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * VcfRemoveGenotypeIfInVcf
 * @author lindenb
 *
 */
public class VcfRemoveGenotypeIfInVcf extends AbstractVCFFilter2 {
	private TabixVcfFileReader tabix=null;
	private boolean removeVariantNoGenotype=false;
	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		long genotypes_reset_count=0;
		long genotypes_empty=0;
		long variant_changed=0;
		long variant_unchanged=0;
		ArrayList<VariantContext> list=new ArrayList<VariantContext>();
		VCFHeader h2=in.getHeader();
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		out.writeHeader(h2);
		while(in.hasNext())
			{
			VariantContext ctx= in.next();
			list.clear();
			Iterator<VariantContext> iter = this.tabix.iterator(
					ctx.getChr(), ctx.getStart(), ctx.getEnd()+1);
			while(iter.hasNext())
				{
				list.add(iter.next());
				}
			
			if(list.isEmpty())
				{
				out.add(ctx);
				++variant_unchanged;
				continue;
				}
			
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			ArrayList<Genotype> genotypes=new ArrayList<Genotype>(ctx.getGenotypes());
			if(genotypes.size()!=h2.getNGenotypeSamples()) throw new RuntimeException();
		
			for(VariantContext snpData: list)
				{
				if(genotypes.isEmpty()) break;
				if(!snpData.getChr().equals(ctx.getChr()) ) continue;
				if(snpData.getStart()!=ctx.getStart() ) continue;
				List<Allele> evsAlts=new ArrayList<Allele>(snpData.getAlternateAlleles());
				if(evsAlts.isEmpty()) continue;
				
				int genotype_index=0;
				while(genotype_index< genotypes.size())
					{
					Genotype g=genotypes.get(genotype_index);

					if(g.isNoCall() || !g.isCalled() || g.isHomRef())
							{
							++genotype_index;
							continue;
							}
					boolean keep=true;
					for(Allele a:g.getAlleles())
						{
						//not same REF anyway
						if(a.isReference() && !a.equals(snpData.getReference())) break;
						
						for(Allele evsAlt:evsAlts)
							{
							//same ALT
							if(a.isNonReference() && a.equals(evsAlt))
								{
								keep=false;
								//warning("Removing Genotype for "+g.getSampleName()+" at "+ctx.getChr()+":"+ctx.getStart());
								genotypes_reset_count++;
								break;
								}
							}
						if(!keep) break;
						}
					if(keep)
						{
						++genotype_index;
						}
					else
						{
						genotypes.remove(genotype_index);
						}
					}	
				}
			//nothing to print
			if(genotypes.isEmpty())
				{
				++genotypes_empty;
				++variant_changed;
				if(removeVariantNoGenotype)
					{
					continue;
					}
				}
			//nothing removed
			if(genotypes.size() == h2.getNGenotypeSamples() )
				{
				out.add(ctx);
				++variant_unchanged;
				continue;
				}
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			++variant_changed;
			}
		info("Number of genotype set to ./. :"+genotypes_reset_count);
		info("Number of variants without any genotype:"+genotypes_empty);
		info("Number of variants changed:"+variant_changed);
		info("Number of variants unchanged:"+variant_unchanged);
		}
	@Override
	public String getProgramDescription() {
		return "Reset Genotypes in VCF if they've been found in another VCF indexed with tabix";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfRemoveGenotypeIfInVcf";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -x (file.vcf.gz) Tabix indexed VCF file. Required");
		out.println(" -r remove variant if there is no called genotype.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String tabixFilePath=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"x:r"))!=-1)
			{
			switch(c)
				{	
				case 'x': tabixFilePath = opt.getOptArg();break;
				case 'r': removeVariantNoGenotype = true;break;
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
		if( tabixFilePath==null)
			{
			error("Undefined VCF tabix  file");
			return -1;
			}
		
		try
			{
			info("Opening "+tabixFilePath);
			this.tabix=new TabixVcfFileReader(tabixFilePath);
			return super.doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.tabix);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfRemoveGenotypeIfInVcf().instanceMainWithExit(args);
	}

}
