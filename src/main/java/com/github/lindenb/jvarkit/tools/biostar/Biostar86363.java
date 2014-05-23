package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class Biostar86363 extends AbstractVCFFilter2
	{
	private static class ChromAndPos
		{
		String chrom;
		int pos;
		ChromAndPos(String chrom,int pos)
			{
			this.chrom=chrom;
			this.pos=pos;
			}
		@Override
		public int hashCode()
			{
			final int prime = 31;
			int result = 1;
			result = prime * result + chrom.hashCode();
			result = prime * result + pos;
			return result;
			}
		@Override
		public boolean equals(Object obj)
			{
			if (this == obj) return true;
			if (obj == null) return false;
			if (!(obj instanceof ChromAndPos)) return false;
			ChromAndPos other = (ChromAndPos) obj;
			if (pos != other.pos) return false;
			return chrom.equals(other.chrom);
			}
		
		}
	private Map<ChromAndPos,Set<String>> pos2sample=new HashMap<Biostar86363.ChromAndPos, Set<String>>();

	private Biostar86363()
		{
		}

	@Override
	public String getProgramDescription()
		{
		return "Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/";
		}
	
	@Override
	public String getProgramName()
		{
		return "Biostar86363";
		}
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar86363";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -h get help (this screen)");
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -G (File) genotypes to reset. Format :CHROM(tab)POS(tab)SAMPLE. REQUIRED.");
		}

	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final List<Allele> empty_g=new ArrayList<Allele>(2);
		empty_g.add(Allele.NO_CALL);
		empty_g.add(Allele.NO_CALL);
		VCFHeader h=in.getHeader();
		final List<String> vcf_samples=h.getSampleNamesInOrder();
		h.addMetaDataLine(new VCFFormatHeaderLine("GR", 1, VCFHeaderLineType.Integer, "(1) = Genotype was reset by "+getProgramName()+":"+getProgramDescription()));
		out.writeHeader(h);
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			ChromAndPos cap=new ChromAndPos(ctx.getChr(), ctx.getStart());
			Set<String> samplesToReset=this.pos2sample.get(cap);
			if(samplesToReset!=null)
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				List<Genotype> genotypes=new ArrayList<Genotype>();
				for(String sample:vcf_samples)
					{
					Genotype g=ctx.getGenotype(sample);
					if(g==null) continue;
					GenotypeBuilder gb=new GenotypeBuilder(g);
					
					if(samplesToReset.contains(sample))
						{
						gb.alleles(empty_g);
						gb.attribute("GR",1);
						}
					else
						{
						gb.attribute("GR",0);
						}
					g=gb.make();
					genotypes.add(g);
					
					}
				vcb.genotypes(genotypes);
				ctx=vcb.make();
				
				}
			out.add(ctx);
			}
		}

	private int loadGenotypes(File f) 
		{
		info("Reading "+f);
		BufferedReader in=null;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			in=IOUtils.openFileForBufferedReading(f);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=tab.split(line);
				if(tokens.length<3)
					{
					error("Bad line in "+line);
					in.close();
					in=null;
					return -1;
					}
				ChromAndPos cap=new ChromAndPos(
					tokens[0],
					Integer.parseInt(tokens[1])
					);
				Set<String> samples=this.pos2sample.get(cap);
				if(samples==null)
					{
					samples=new HashSet<String>();
					this.pos2sample.put(cap,samples);
					}
				samples.add(tokens[2]);
				}
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
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args, "hvL:G:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(opt.getOptArg()));break;
				case 'G': if(loadGenotypes(new File(opt.getOptArg()))!=0) return -1;break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+opt.getOptOpt());return -1;
				}
			}
		return doWork(opt.getOptInd(), args);
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar86363().instanceMainWithExit(args);
		}

	}
