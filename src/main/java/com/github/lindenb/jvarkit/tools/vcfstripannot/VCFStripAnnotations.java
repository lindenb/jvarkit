package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;


import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VCFStripAnnotations extends AbstractVCFFilter2
	{
	private Set<String> KEY=new HashSet<String>();
	private Set<String> FORMAT=new HashSet<String>();
	private Set<String> FILTER=new HashSet<String>();
	
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VCFStripAnnotations";
		}
	
	@Override
	protected String getProgramCommandLine()
		{
		return " Removes one or more field from the INFO/FORMAT column of a VCF.";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -k (key) remove this INFO attribute. '*'= all keys");
		out.println(" -f (format) remove this FORMAT attribute. '*'= all keys BUT GT/DP/AD/GQ/PL");
		out.println(" -F (filter) remove this FILTER. '*'= all keys.");
		super.printOptions(out);
		}
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
			{
			VCFHeader header=r.getHeader();
			boolean remove_all_info=this.KEY.contains("*");
			boolean remove_all_format=this.FORMAT.contains("*");
			boolean remove_all_filters=this.FILTER.contains("*");
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			
			for(Iterator<VCFInfoHeaderLine> h=h2.getInfoHeaderLines().iterator();
					h.hasNext();)
				{
				VCFInfoHeaderLine vih=h.next();
				if(this.KEY.contains(vih.getID()))
					h.remove();
				}
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			w.writeHeader(h2);
			
			while(r.hasNext())
				{
				VariantContext ctx=r.next();
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				/* INFO */
				if(remove_all_info)
					{
					for(String k2: ctx.getAttributes().keySet())
						{
						b.rmAttribute(k2);
						}
					}
				else if(!KEY.isEmpty())
					{
					for(String key:KEY) b.rmAttribute(key);
					}
				
				/* formats */
				if(remove_all_format)
					{
					List<Genotype> genotypes=new ArrayList<Genotype>();
					for(Genotype g:ctx.getGenotypes())
						{
						GenotypeBuilder gb=new GenotypeBuilder(g);
						gb.attributes(new HashMap<String, Object>());
						genotypes.add(gb.make());
						}
					b.genotypes(genotypes);
					}
				else if(! this.FORMAT.isEmpty())
					{
					List<Genotype> genotypes=new ArrayList<Genotype>();
					for(Genotype g:ctx.getGenotypes())
						{
						GenotypeBuilder gb=new GenotypeBuilder(g);
						Map<String, Object> map=new HashMap<String, Object>();
						for(String key: g.getExtendedAttributes().keySet())
							{
							if(this.FORMAT.contains(key)) continue;
							map.put(key, g.getExtendedAttribute(key));
							}
						gb.attributes(map);
						genotypes.add(gb.make());
						}
					b.genotypes(genotypes);
					}
				
				/* filters */
				if(remove_all_filters)
					{
					b.unfiltered();
					}
				else if(! this.FILTER.isEmpty())
					{
					b.unfiltered();
					for(String key:ctx.getFilters())
						{
						if(this.FILTER.contains(key)) continue;
						if(key.equals("PASS")) continue;
						b.filter(key);
						}
					}
				w.add(b.make());
				}		
			}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "k:F:f:"))!=-1)
			{
			switch(c)
				{
				case 'k': this.KEY.add(opt.getOptArg()); break;
				case 'f': this.FORMAT.add(opt.getOptArg()); break;
				case 'F': this.FILTER.add(opt.getOptArg()); break;
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

	
	public static void main(String[] args) throws IOException
		{
		new VCFStripAnnotations().instanceMainWithExit(args);
		}
	}
