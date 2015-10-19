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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FindAVariation extends AbstractFindAVariation
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FindAVariation.class);

	
	private static class Mutation
		{
		String chrom;
		int pos;
		Mutation(String chrom,int pos)
			{
			this.chrom=chrom;
			this.pos=pos;
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + chrom.hashCode();
			result = prime * result + pos;
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)return true;
			Mutation other = (Mutation) obj;
			if (pos != other.pos) return false;
			 if (!chrom.equals(other.chrom))
				return false;
			
			return true;
		}
		
		@Override
		public String toString() {
			return chrom+":"+pos;
			}
		
		}
	
	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFindAVariation.AbstractFindAVariationCommand
	 	{		
		private Set<Mutation> mutations=new HashSet<Mutation>();
		private PrintWriter out=null;
	
   
    
    private void reportPos(File f,VariantContext ctx)
		{
		out.print(f);
		out.print('\t');
		out.print(ctx.getContig());
		out.print('\t');
		out.print(ctx.getStart());
		out.print('\t');
		out.print(ctx.getEnd());
		out.print('\t');
		out.print(ctx.hasID()?ctx.getID():".");
		out.print('\t');
		}	

    
    private void report(File f,VariantContext ctx)
    	{
    	GenotypesContext genotypes=ctx.getGenotypes();
    	if(genotypes==null || genotypes.isEmpty())
    		{
    		reportPos(f,ctx);
    		out.println();
    		}
    	else
    		{
    		for(int i=0;i< genotypes.size();++i)
    			{
    			Genotype g=genotypes.get(i);
    			reportPos(f,ctx);
    			out.print('\t');
    			out.print(g.getSampleName());
    			out.print('\t');
    			out.print(g.getType());
    			out.print('\t');
    			List<Allele> alleles=g.getAlleles();
    			for(int na=0;na<alleles.size();++na)
    				{
    				if(na>0) out.print(" ");
    				out.print(alleles.get(na).getDisplayString());
    				}
    			out.println();
    			}
    		}
    	}	
    
    private Set<Mutation> convertFromVcfHeader(File f,VCFHeader h)
    	{
    	Set<Mutation> copy=new HashSet<Mutation>(this.mutations.size());
    	for(Mutation m:this.mutations)
    		{
    		String s=VCFUtils.findChromNameEquivalent(m.chrom,h);
    		if(s==null)
    			{
    			LOG.warn("Cannot convert chrom "+s+" in "+f);
    			continue;
    			}
    		copy.add(new Mutation(s, m.pos));
    		}
    	return copy;
    	}

    private void scan(BufferedReader in) throws IOException
    	{
    	
    	String line;
    	while((line=in.readLine())!=null)
    			{
    			if(line.isEmpty() || line.startsWith("#")) continue;
    			File f=new File(line);
    			if(!f.isFile()) continue;
    			if(!f.canRead()) continue;
    			if(!VCFUtils.isVcfFile(f)) continue;
    			LOG.info(f);
    			VcfIterator iter=null;
    			
	    			if(VCFUtils.isTabixVcfFile(f))
	    				{
	    				TabixVcfFileReader r=null;
		    			try
							{
							r=new TabixVcfFileReader(f.getPath());
							for(Mutation m:convertFromVcfHeader(f,r.getHeader()))
								{
								Iterator<VariantContext> iter2=r.iterator(
										m.chrom, m.pos, m.pos);
								while(iter2.hasNext())
									{
									report(f,iter2.next());
									}
								CloserUtil.close(iter2);
								}
							}
		    			catch(htsjdk.tribble.TribbleException.InvalidHeader err)
		    				{
		    				LOG.warn(f+"\t"+err.getMessage());
		    				}
						catch(Exception err)
							{
							LOG.error(err);
							}
						finally
							{
							CloserUtil.close(r);
							}    				
	    				}
	    			else
	    				{
	    				try
	    					{
	    					iter=VCFUtils.createVcfIteratorFromFile(f);
	    					Set<Mutation> mutlist=convertFromVcfHeader(f,iter.getHeader());
	    					while(iter.hasNext())
	    						{
	    						VariantContext ctx=iter.next();
	    						Mutation m=new Mutation(ctx.getContig(), ctx.getStart());
	    						if(mutlist.contains(m))
	    							{
	    							report(f,ctx);
	    							}
	    						}
	    					
	    					}
	    				catch(htsjdk.tribble.TribbleException.InvalidHeader err)
		    				{
		    				LOG.warn(f+"\t"+err.getMessage());
		    				}
	    				catch(Exception err)
	    					{
	    					LOG.error(err);
	    					}
	    				finally
	    					{
	    					CloserUtil.close(iter);
	    					}
	    				}
	    			
    			}
    	}
    @Override
    public Collection<Throwable> call() throws Exception
    	{
    	if(super.positionStrSet.isEmpty())
			{
			return wrapException("position not defined");
			}
		else for(String s: super.positionStrSet)
			{
			int colon=s.indexOf(':');
			if(colon==-1 || colon+1==s.length())
				{
				return wrapException("Bad chrom:pos "+s);
				}
			
			String chrom=s.substring(0,colon).trim();
			if(chrom.isEmpty())
				{
				return wrapException("Bad chrom:pos "+s);
				}
			Mutation m=new Mutation(chrom, Integer.parseInt(s.substring(colon+1)));
			LOG.info("mutation =  "+m);
			this.mutations.add(m);
			}
		if(this.mutations.isEmpty())
			{
			return wrapException("position not defined");
			}
		final List<String> args = getInputFiles();
		this.out = null;
		BufferedReader r = null;
		try
			{
			this.out =openFileOrStdoutAsPrintWriter();  
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				 r=new BufferedReader(new InputStreamReader(stdin()));
				this.scan(r);
				r.close();
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("Reading from "+filename);
					r=IOUtils.openURIForBufferedReading(filename);
					this.scan(r);
					r.close();
					}
				}
			this.out.flush();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(this.out);
			CloserUtil.close(r);
			}
		}
	 	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindAVariation().instanceMainWithExit(args);

	}

}
