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
package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
 * 
 * VCFFixIndels
 *
 */
public class VCFFixIndels extends AbstractVCFFilter3
	{
	public VCFFixIndels()
		{
		}
	
	@Override
	public String getProgramDescription() {
		return " Fix samtools indels (for @SolenaLS).";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VCFFixIndels";
		}
	
	
	@Override
	protected void doWork(
			String source,
			VcfIterator r,
			VariantContextWriter w
			)
			throws IOException
		{
		long nChanged=0L;
		final String TAG="INDELFIXED";
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,VCFHeaderLineType.String,"Fix Indels for @SolenaLS (position|alleles...)"));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		w.writeHeader(h2);
		SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
		
		while(r.hasNext())
			{
			boolean somethingChanged=false;
			VariantContext ctx=progress.watch(r.next());
			/* map old allele to new allele */
			Map<Allele, Allele> original2modified=new HashMap<Allele, Allele>();
			/* did we found a strange allele (symbolic, etc ) */
			boolean strange_allele_flag=false;
			for(Allele a:ctx.getAlleles())
				{
				original2modified.put(a, a);
				if(a.isSymbolic() || a.isNoCall())
					{
					strange_allele_flag=true;
					break;
					}
				}
			
			if(strange_allele_flag ||
				original2modified.size()<2 /* at least 2 alleles: REF+ ALT */
				)
				{
				incrVariantCount();
				w.add(ctx);
				continue;
				}
			
			/* record chromStart if we need to shift to the right */
			int chromStart= ctx.getStart();
			
			/* trim right then left  */
			for(int side=0;side<2;++side)
				{
				boolean done=false;
				while(!done)
					{
					boolean ok_side=true;
					/* common nucleotide at end/start */
					Character targetChar = null;
					done=true;
					//scan side
					Set<Allele> keys = new HashSet<>(original2modified.keySet());  
					for(Allele k:keys)
						{
						Allele newAllele = original2modified.get(k);
						if(newAllele.isSymbolic()) 
							{
							ok_side = false;
							break;
							}
						String baseString = newAllele.getBaseString().trim().toUpperCase();
						if(baseString.length()<2)
							{
							ok_side = false;
							break;
							}
						/* first or last char or all sequences
						 * side==0 : right
						 * side==1 : left
						 * */
						Character baseChar=
							(side==0?
							baseString.charAt(baseString.length()-1):
							baseString.charAt(0)
							);
						if(targetChar==null)
							{
							targetChar = baseChar;
							}
						else if(!targetChar.equals(baseChar))
							{
							/* doesn't end with same nucleotide */
							ok_side = false;
							break;
							}
						}
					/* ok we can shift all alleles */
					if(ok_side && targetChar!=null)
						{
						done=false;
						somethingChanged=true;
						for(Allele k:keys)
							{
							Allele newAllele = original2modified.get(k);
							String baseString = newAllele.getBaseString().trim().toUpperCase();
							if(side==0)//remove last nucleotide
								{
								newAllele = Allele.create(
									baseString.substring(0,baseString.length()-1),
									newAllele.isReference());
								}
							else
								{
								newAllele = Allele.create(
										baseString.substring(1),
										newAllele.isReference());
								
								}
							original2modified.put(k, newAllele);
							}
						if(side==1) chromStart++;
						}
					}/* end of while done */
				
				}/* end side */
			
			if(!somethingChanged)
				{
				w.add(ctx);
				incrVariantCount();
				continue;
				}
			
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			b.start(chromStart);
			Allele newRef=original2modified.get(ctx.getReference());
			b.stop(chromStart+newRef.getBaseString().length()-1);
			b.alleles(original2modified.values());
			List<Genotype> genotypes=new ArrayList<>();
			for(String sample:header.getSampleNamesInOrder())
				{
				Genotype g = ctx.getGenotype(sample);
				if(g.isNoCall()) 
					{
					genotypes.add(g);
					continue;
					}
				GenotypeBuilder gb=new GenotypeBuilder(g);
				List<Allele> aL=new ArrayList<>();
				for(Allele a:g.getAlleles())
					{
					aL.add(original2modified.get(a));
					}
				gb.alleles(aL);
				genotypes.add(gb.make());
				}
			
			StringBuilder tagContent=new StringBuilder();
			tagContent.append(String.valueOf(ctx.getStart()));
			for(Allele a:ctx.getAlleles())
				{
				tagContent.append("|");
				tagContent.append(a.toString());
				}
			
			
			b.attribute(TAG,tagContent.toString());
			b.genotypes(genotypes);
			
			incrVariantCount();
			w.add( b.make());
			
			if(this.checkOutputError()) break;
			}
		progress.finish();
		info("indels changed:"+nChanged);
		}
	

	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-o (filename) file out. default:stdout");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(opt.getOptArg());break;
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
		return super.mainWork(opt.getOptInd(), args);
		}
	
	public static void main(String[] args) throws IOException
		{
		new VCFFixIndels().instanceMainWithExit(args);
		}
}
