package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class Biostar86363 extends AbstractBiostar86363
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar86363.class);
	
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
	@Override
	public Command createCommand() {
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar86363.AbstractBiostar86363Command
		{
	private Map<ChromAndPos,Set<String>> pos2sample=new HashMap<Biostar86363.ChromAndPos, Set<String>>();


	
	
	private void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final List<Allele> empty_g=new ArrayList<Allele>(2);
		empty_g.add(Allele.NO_CALL);
		empty_g.add(Allele.NO_CALL);
		VCFHeader h=in.getHeader();
		final List<String> vcf_samples=h.getSampleNamesInOrder();
		h.addMetaDataLine(new VCFFormatHeaderLine("GR", 1, VCFHeaderLineType.Integer, "(1) = Genotype was reset by "+getName()+":"+getFactory().getDescription()));
		out.writeHeader(h);
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			ChromAndPos cap=new ChromAndPos(ctx.getContig(), ctx.getStart());
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

	private int loadGenotypes(String uri) throws IOException
		{
		LOG.info("Reading "+uri);
		BufferedReader in=null;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			in=IOUtils.openURIForBufferedReading(uri);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=tab.split(line);
				if(tokens.length<3)
					{
					in.close();
					in=null;
					throw new IOException("Bad line in "+line);
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
		finally
			{
			CloserUtil.close(in);
			}	
		}
	
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			if( genotypeURI == null)
				{
				return wrapException("undefined list of genotypes");
				}
			
			VcfIterator in=null;
			VariantContextWriter out=null;
			try {
				loadGenotypes(genotypeURI);
				in  = super.openVcfIterator(inputName);
				out = super.openVariantContextWriter();
				doWork(in, out);
				return Collections.emptyList();
			} catch (Exception e) {
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(in);
				CloserUtil.close(out);
				}
			}
	
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar86363().instanceMainWithExit(args);
		}

	}
