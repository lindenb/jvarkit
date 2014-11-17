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

import java.io.PrintWriter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class AlleleFrequencyCalculator extends AbstractCommandLineProgram
	{
	private AlleleFrequencyCalculator()
		{
		
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/AlleleFrequencyCalculator";
		}
	
	@Override
	public String getProgramDescription() {
		return "Allele Frequency Calculator";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
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
		
		VcfIterator in=null;
		PrintWriter out=new PrintWriter(System.out);
		try
			{
			if(opt.getOptInd()==args.length)
				{
				in=VCFUtils.createVcfIteratorStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				in=VCFUtils.createVcfIterator(args[opt.getOptInd()]);
				}
			else
				{
				error(super.getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			out.println("CHR\tPOS\tID\tREF\tALT\tTOTAL_CNT\tALT_CNT\tFRQ");
			while(in.hasNext() && !out.checkError())
				{
				
				VariantContext ctx=in.next();
				Allele ref=ctx.getReference();
				if(ref==null) continue;
				if(ctx.getNSamples()==0 || ctx.getAlternateAlleles().isEmpty()) continue;
				Allele alt=ctx.getAltAlleleWithHighestAlleleCount();
				if(alt==null) continue;
				
				GenotypesContext genotypes=ctx.getGenotypes();
				if(genotypes==null) continue;
				int total_ctn=0;
				int alt_ctn=0;
				for(int i=0;i< genotypes.size();++i)
					{
					Genotype g=genotypes.get(i);
					for(Allele allele: g.getAlleles())
						{
						if(allele.equals(ref))
							{
							total_ctn++;
							}
						else if (allele.equals(alt))
							{
							total_ctn++;
							alt_ctn++;
							}
						}
					
					}
				
				
				out.print(ctx.getChr());
				out.print("\t");
				out.print(ctx.getStart());
				out.print("\t");
				out.print(ctx.hasID()?ctx.getID():".");
				out.print("\t");
				out.print(ref.getBaseString());
				out.print("\t");
				out.print(alt.getBaseString());
				out.print("\t");
				out.print(total_ctn);
				out.print("\t");
				out.print(alt_ctn);
				out.print("\t");
				out.print(alt_ctn/(float)total_ctn);
				out.println();
				}
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	public static void main(String[] args) {
		new AlleleFrequencyCalculator().instanceMainWithExit(args);

	}

}
