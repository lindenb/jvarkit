/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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


*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

@Program(name="vcfgroupbypop",
	description="Group VCF data by population, creates a VCF  where each 'SAMPLE' is a population",
	creationDate="20190319",
	modificationDate="20210528",
	keywords= {"vcf","pedigree","population"}
	)
public class VcfGroupByPopulation extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfGroupByPopulation.class).make();
	@Parameter(names={"-p","--mapping"},description="mapping file: each line is (SAMPLE)\\t(POP)\\n",required=true)
	private Path mappingFile=null;
	@Parameter(names={"-m","--min-fisher"},description="min inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double min_fisher = 0.0 ;
	@Parameter(names={"-M","--max-fisher"},description="max inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double max_fisher = 1.0 ;

	private final Map<String,String> sample2population=new TreeMap<>();
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private static class GCount
		{
		int R=0;
		int A=0;
		int uncalled=0;
		int dp=-1;
		int have_alt = 0;
		int miss_alt = 0;
		private void watch(final Genotype g)
			{
			if(!g.isAvailable() ||
			   !g.isCalled() ||
			    g.isNoCall())
				{
				++uncalled;
				return;
				}
			if(g.isHomRef())
				{
				R+=2;
				++miss_alt;
				}
			else if(g.isHetNonRef())
				{
				A+=2;
				++have_alt;
				}
			else if(g.isHet())
				{
				R++;
				A++;
				++have_alt;
				}
			else if(g.isHomVar())
				{
				A+=2;
				++have_alt;
				}
			if(g.hasDP())
				{
				if(this.dp==-1) dp=0;
				this.dp+=g.getDP();
				}
			}
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator vcfIn, VariantContextWriter out) {
			try {
			final VCFHeader header= vcfIn.getHeader();
			final Set<String> samplesInVcf=new HashSet<>( header.getSampleNamesInOrder());
			
			this.sample2population.keySet().retainAll(samplesInVcf);
			
			final Map<String,Set<String>> population2samples=new TreeMap<>();
			for(String sample:this.sample2population.keySet())
				{
				final String pop= this.sample2population.get(sample);
				Set<String> samples= population2samples.get(pop);
				if(samples==null)
					{
					samples=new HashSet<>();
					population2samples.put(pop,samples);
					}
				samples.add(sample);
				}
			final List<String> poplist = new ArrayList<>(population2samples.keySet());
			
			for(String sample: header.getSampleNamesInOrder())
				{
				if(!this.sample2population.containsKey(sample))
					{
					throw new IOException("Sample "+sample+" not affected to a population");
					}
				}
			
			final Set<VCFHeaderLine> metaData=new LinkedHashSet<>();
	
			for(int i=0;i< poplist.size();i++) {
					for(int j=i+1;j< poplist.size();j++) {
					metaData.add( new VCFInfoHeaderLine(
							poplist.get(i)+"_"+poplist.get(j)+"_FISHER",1,VCFHeaderLineType.Float,"Fisher Test"
							));
					}
				}
			
			/* FORMAT */
			metaData.add( new VCFFormatHeaderLine(
					"NS",1,VCFHeaderLineType.Integer,"Total Number of Samples"
					));
	
			metaData.add( new VCFFormatHeaderLine(
					"R",1,VCFHeaderLineType.Integer,"Number of alleles REF (hom:=2,het:=1)"
					));
			metaData.add( new VCFFormatHeaderLine(
					"A",1,VCFHeaderLineType.Integer,"Number of alleles ALT (hom:=2,het:=1)"
					));
			metaData.add( new VCFFormatHeaderLine(
					"UNC",1,VCFHeaderLineType.Integer,"Number of NON-called samples"
					));
			metaData.add( new VCFFormatHeaderLine(
					"F",1,VCFHeaderLineType.Float,"Allele Frequency A/(R+A)"
					));
	
			metaData.add(new VCFFormatHeaderLine(
					"DP",
					1,
					VCFHeaderLineType.Integer,
					"Depth"));
	
			
			/* INFO */
			metaData.add( new VCFInfoHeaderLine(
					"NS",1,VCFHeaderLineType.Integer,"Total Number of Samples"
					));
	
			metaData.add( new VCFInfoHeaderLine(
					"R",1,VCFHeaderLineType.Integer,"Number of alleles REF (hom:=2,het:=1)"
					));
			metaData.add( new VCFInfoHeaderLine(
					"A",1,VCFHeaderLineType.Integer,"Number of alleles ALT (hom:=2,het:=1)"
					));
			metaData.add( new VCFInfoHeaderLine(
					"UNC",1,VCFHeaderLineType.Integer,"Number of NON-called samples"
					));
			metaData.add( new VCFInfoHeaderLine(
					"F",1,VCFHeaderLineType.Float,"Allele Frequency A/(R+A)"
					));
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.DEPTH_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.GENOTYPE_KEY);
	
			
			final VCFHeader h2=new VCFHeader(
					metaData,
					population2samples.keySet()
					);
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null) h2.setSequenceDictionary(dict);

			
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
					
			while(vcfIn.hasNext())
				{
				final VariantContext ctx= vcfIn.next();
				final VariantContextBuilder vcb=new VariantContextBuilder(
						ctx.getSource(),
						ctx.getContig(),
						ctx.getStart(),
						ctx.getEnd(),
						ctx.getAlleles()
						);
				if(ctx.hasID()) vcb.id(ctx.getID());
				GCount count_ctx = new GCount();
				final Map<String,GCount> pop2gcount = new HashMap<>();
				final List<Genotype> genotypes= new ArrayList<>(population2samples.size());
				for(String pop:population2samples.keySet())
					{
					final GCount count=new GCount();
					final Set<String> samples = population2samples.get(pop);
					pop2gcount.put(pop, count);
					for(String sample: samples)
						{
						final Genotype g= ctx.getGenotype(sample);
						count.watch(g);
						}
					
					final GenotypeBuilder gb=new GenotypeBuilder(pop);
					
					gb.attribute("NS", samples.size());
					gb.attribute("R", count.R);
					gb.attribute("A", count.A);
					gb.attribute("UNC", count.uncalled);
					if(count.R+count.A>0)
						{
						gb.attribute("F",
								(float)count.A/(float)(count.R+count.A)
								);
						}
					if(count.dp>=0)
						{
						gb.attribute("DP", count.dp);
						if(count_ctx.dp==-1) count_ctx.dp=0;
						}
					
					genotypes.add(gb.make());
					
					count_ctx.R += count.R;
					count_ctx.A += count.A;
					count_ctx.uncalled += count.uncalled;
					count_ctx.dp += count.dp;
	
					}
				vcb.attribute("R", count_ctx.R);
				vcb.attribute("A", count_ctx.A);
				vcb.attribute("UNC", count_ctx.uncalled);
				if(count_ctx.R+count_ctx.A>0)
					{
					vcb.attribute("F",
							(float)count_ctx.A/(float)(count_ctx.R+count_ctx.A)
							);
					}
				if(count_ctx.dp>=0)
					{
					vcb.attribute(VCFConstants.DEPTH_KEY, count_ctx.dp);
					}
				vcb.attribute("NS", this.sample2population.keySet().size());
				
				
				for(int i=0;i< poplist.size();i++) {
					final GCount gi = pop2gcount.get(poplist.get(i));
					for(int j=i+1;j< poplist.size();j++) {
							final GCount gj = pop2gcount.get(poplist.get(j));
							final FisherExactTest fisherTest = FisherExactTest.compute(
									gi.have_alt,gi.miss_alt,
									gj.have_alt,gj.miss_alt
									);
							double fisher = fisherTest.getAsDouble();
							if(fisher> this.max_fisher) continue;
							if(fisher< this.min_fisher) continue;
							vcb.attribute(poplist.get(i)+"_"+poplist.get(j)+"_FISHER",fisher);
					}
				}

				
				vcb.genotypes(genotypes);
				out.add(vcb.make());
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected int beforeVcf()
		{
		try(BufferedReader br = IOUtils.openPathForBufferedReading(this.mappingFile)) {
			String line;
			while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
					int space=line.indexOf('\t');
					if(space<=0) throw new IOException("tab missing in "+line);
					String sample = line.substring(0,space);
					String pop= line.substring(space+1);
					if(StringUtils.isBlank(sample))  throw new IOException("empty sample in "+line);
					if(StringUtils.isBlank(pop))  throw new IOException("empty sample in "+line);
					String prevpo= this.sample2population.get(sample);
					if(prevpo!=null && !prevpo.equals(pop))
						throw new IOException("two pop declared for "+sample);
					this.sample2population.put(sample, pop);
					}
			}
		catch(final IOException err ) {
			LOG.error(err);
			return -1;
			}
		
		return super.beforeVcf();
		}
	
	
	public static void main(final String[] args)
		{
		new VcfGroupByPopulation().instanceMainWithExit(args);
		}
}
