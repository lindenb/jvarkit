package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

@Program(name="biostar86363",description="Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/")

public class Biostar86363 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar86363.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	
	
	private Map<ContigPosRef,Set<String>> pos2sample=new HashMap<ContigPosRef, Set<String>>();

	private Biostar86363()
		{
		}

	@Parameter(names="-G",description="genotypes to reset. Format :CHROM(tab)POS(tab)ref(tab)SAMPLE. REQUIRED.",required=true)
	private File genotypeFile=null;
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) 
		{
		final List<Allele> empty_g=new ArrayList<Allele>(2);
		empty_g.add(Allele.NO_CALL);
		empty_g.add(Allele.NO_CALL);
		VCFHeader h=in.getHeader();
		final List<String> vcf_samples=h.getSampleNamesInOrder();
		h.addMetaDataLine(new VCFFormatHeaderLine("GR", 1, VCFHeaderLineType.Integer, "(1) = Genotype was reset by "+getProgramName()));
		out.writeHeader(h);
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			ContigPosRef cap=new ContigPosRef(ctx);
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
				ctx=VCFUtils.recalculateAttributes(vcb.make());
				}
			out.add(ctx);
			}
		return 0;
		}

	
	@Override
	public int doWork(List<String> args) {
		if(genotypeFile==null)
			{
			LOG.error("undefined genotype file");
			return -1;
			}
		BufferedReader in=null;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			in=IOUtils.openFileForBufferedReading(this.genotypeFile);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=tab.split(line);
				if(tokens.length<4)
					{
					LOG.error("Bad line in "+line);
					in.close();
					in=null;
					return -1;
					}
				ContigPosRef cap=new ContigPosRef(
					tokens[0],
					Integer.parseInt(tokens[1]),
					Allele.create(tokens[2],true)
					);
				Set<String> samples=this.pos2sample.get(cap);
				if(samples==null)
					{
					samples=new HashSet<String>();
					this.pos2sample.put(cap,samples);
					}
				samples.add(tokens[3]);
				}
			in.close();
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(in);
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
