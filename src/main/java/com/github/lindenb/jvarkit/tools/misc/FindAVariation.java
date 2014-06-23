package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
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
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FindAVariation extends AbstractCommandLineProgram
	{
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
	private Set<Mutation> mutations=new HashSet<Mutation>();
	private PrintWriter out=null;
	
    private FindAVariation()
    	{
    	
    	}		
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/FindAVariation";
    	}
    
    @Override
    public String getProgramDescription() {
    	return "Find a specific mutation in a list of VCF files.";
    	}
    
    private void reportPos(File f,VariantContext ctx)
		{
		out.print(f);
		out.print('\t');
		out.print(ctx.getChr());
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
    			warning("Cannot convert chrom "+s+" in "+f);
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
    			info(f);
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
		    				warning(f+"\t"+err.getMessage());
		    				}
						catch(Exception err)
							{
							error(err);
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
	    						Mutation m=new Mutation(ctx.getChr(), ctx.getStart());
	    						if(mutlist.contains(m))
	    							{
	    							report(f,ctx);
	    							}
	    						}
	    					
	    					}
	    				catch(htsjdk.tribble.TribbleException.InvalidHeader err)
		    				{
		    				warning(f+"\t"+err.getMessage());
		    				}
	    				catch(Exception err)
	    					{
	    					error(err);
	    					}
	    				finally
	    					{
	    					CloserUtil.close(iter);
	    					}
	    				}
	    			
    			}
    	}
    
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -p chrom:pos . Add this chrom/position.");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"p:"))!=-1)
			{
			switch(c)
				{
				case 'p':
					{
					String s=opt.getOptArg();
					int colon=s.indexOf(':');
					if(colon==-1 || colon+1==s.length())
						{
						error("Bad chrom:pos "+s);
						return -1;
						}
					
					String chrom=s.substring(0,colon).trim();
					if(chrom.isEmpty())
						{
						error("Bad chrom:pos "+s);
						return -1;
						}
					Mutation m=new Mutation(chrom, Integer.parseInt(s.substring(colon+1)));
					info("adding "+m);
					this.mutations.add(m);
					break;
					}
				default:
					{
					switch(super.handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS:return 0;
						case OK:break;
						}
					}
				}
			}
		try
			{
			this.out=new PrintWriter(System.out);
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				scan(new BufferedReader(new InputStreamReader(System.in)));
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					BufferedReader r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();
					}
				}
			this.out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{

			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindAVariation().instanceMainWithExit(args);

	}

}
