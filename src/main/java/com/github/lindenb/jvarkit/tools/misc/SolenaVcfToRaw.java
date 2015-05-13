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

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.List;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;



public class SolenaVcfToRaw extends AbstractKnimeApplication
	{
	
	public SolenaVcfToRaw()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Solena: vcf to (chrom/pos/ref/alt/individus(G:0/1/2/-9)";
		}
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"SolenaVcfToRaw";
		}
	
	@Override
	public int initializeKnime() {
		return super.initializeKnime();
		}

	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		PrintWriter pw=null;
		VcfIterator in=null;
		try
			{
			if(args.isEmpty())
				{
				info("reading from stdin.");
				in = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				info("reading from "+filename);
				in = VCFUtils.createVcfIterator(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			if(getOutputFile()==null)
				{
				pw = new PrintWriter(System.out);
				}
			else
				{
				pw = new PrintWriter(getOutputFile());
				}
			List<String> samples= in.getHeader().getSampleNamesInOrder();
			pw.print("CHROM\tPOS\tREF\tALT");
			for(String sample:samples)
				{
				pw.print("\t");
				pw.print(sample);
				}
			pw.println();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
			while(in.hasNext())
				{
				VariantContext ctx= progress.watch(in.next());
				if(ctx.getAlternateAlleles().size()!=1)
					{
					info("count(ALT)!=1 in "+ctx.getChr()+":"+ctx.getStart());
					continue;
					}
				pw.print(ctx.getChr());
				pw.print("\t");
				pw.print(ctx.getStart());
				pw.print("\t");
				pw.print(ctx.getReference().getDisplayString());
				pw.print("\t");
				pw.print(ctx.getAlternateAlleles().get(0).getDisplayString());
				for(String sample:samples)
					{
					Genotype g=ctx.getGenotype(sample);
					pw.print("\t");
					if(g.isHomRef())
						{
						pw.print("0");
						}
					else if(g.isHomVar())
						{
						pw.print("2");
						}
					else if(g.isHet())
						{
						pw.print("1");
						}
					else
						{
						pw.print("-9");
						}
					}
				pw.println();
				if(checkOutputError()) break;
				}
			progress.finish();
			
			
			pw.flush();
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
			CloserUtil.close(pw);
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (file) output file (default stdout)");
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
				case 'o': setOutputFile(opt.getOptArg()); break;
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
		new SolenaVcfToRaw().instanceMainWithExit(args);
		}
	}
