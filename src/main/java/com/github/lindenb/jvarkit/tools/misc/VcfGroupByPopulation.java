/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

@Program(name="vcfgroupbypop",
description="Group VCF data by population, creates a VCF  where each 'SAMPLE' is a population")
public class VcfGroupByPopulation extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGroupByPopulation.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	@Parameter(names="-p",description="mapping file: each line is (SAMPLE)\\t(POP)\\n",required=true)
	private File mappingFile=null;


	private Map<String,String> sample2population=new HashMap<>();
	
	private static class GCount
		{
		int R=0;
		int A=0;
		int uncalled=0;
		int dp=-1;
		private void watch(Genotype g)
			{
			if(!g.isAvailable() ||
			   !g.isCalled() ||
			    g.isNonInformative() ||
			    g.isNoCall())
				{
				++uncalled;
				return;
				}
			if(g.isHomRef())
				{
				R+=2;
				}
			else if(g.isHetNonRef())
				{
				A+=2;
				}
			else if(g.isHet())
				{
				R++;
				A++;
				}
			else if(g.isHomVar())
				{
				A+=2;
				}
			if(g.hasDP())
				{
				if(this.dp==-1) dp=0;
				this.dp+=g.getDP();
				}
			}
		}
	
	public VcfGroupByPopulation()
		{
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator vcfIn, VariantContextWriter out) {
		try {
			
			Reader r=IOUtils.openFileForBufferedReading(this.mappingFile);
			parsePopulationMapping(r);
			r.close();
			
		VCFHeader header= vcfIn.getHeader();
		Set<String> samplesInVcf=new HashSet<>( header.getSampleNamesInOrder());
		
		this.sample2population.keySet().retainAll(samplesInVcf);
		
		Map<String,Set<String>> population2samples=new HashMap<>();
		for(String sample:this.sample2population.keySet())
			{
			String pop= this.sample2population.get(sample);
			Set<String> samples= population2samples.get(pop);
			if(samples==null)
				{
				samples=new HashSet<>();
				population2samples.put(pop,samples);
				}
			samples.add(sample);
			}
		
		for(String sample: header.getSampleNamesInOrder())
			{
			if(!this.sample2population.containsKey(sample))
				{
				throw new IOException("Sample "+sample+" not affected to a population");
				}
			}
		
		Set<VCFHeaderLine> metaData=new LinkedHashSet<>();
		for(VCFHeaderLine vhl: header.getMetaDataInInputOrder())
			{
			if(!(vhl instanceof VCFContigHeaderLine)) continue;
			metaData.add(vhl);
			}
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		metaData.addAll(JVarkitVersion.getInstance().getMetaData(getClass().getSimpleName()));

		
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
		metaData.add(new VCFInfoHeaderLine(
				"DP",
				1,
				VCFHeaderLineType.Integer,
				"Depth"));

		
		
		metaData.add(new VCFFormatHeaderLine(
				"GT",
				1,
				VCFHeaderLineType.String,
				"Genotype"));
		VCFHeader h2=new VCFHeader(
				metaData,
				population2samples.keySet()
				);
		
		out.writeHeader(h2);
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(vcfIn.getHeader());
		
		while(vcfIn.hasNext())
			{
			VariantContext ctx=progress.watch(vcfIn.next());
			VariantContextBuilder vcb=new VariantContextBuilder(
					ctx.getSource(),
					ctx.getContig(),
					ctx.getStart(),
					ctx.getEnd(),
					ctx.getAlleles()
					);
			if(ctx.hasID()) vcb.id(ctx.getID());
			GCount count_ctx = new GCount();
			
			List<Genotype> genotypes= new ArrayList<>(population2samples.size());
			for(String pop:population2samples.keySet())
				{
				GCount count=new GCount();
				Set<String> samples = population2samples.get(pop);
				for(String sample: samples)
					{
					Genotype g= ctx.getGenotype(sample);
					count.watch(g);
					}
				
				GenotypeBuilder gb=new GenotypeBuilder(pop);
				
				
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
				vcb.attribute("DP", count_ctx.dp);
				}
			vcb.attribute("NS", this.sample2population.keySet().size());
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			}
		progress.finish();
		return 0;
		}
	catch(Exception err) {
		LOG.error(err);
		return -1;
		}
	}
	
	
	
	public void parsePopulationMapping(Reader in) throws IOException
		{
		BufferedReader r=new BufferedReader(in);
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			int space=line.indexOf('\t');
			if(space<=0) throw new IOException("tab missing in "+line);
			String sample = line.substring(0,space);
			String pop= line.substring(space+1);
			if(sample.trim().isEmpty())  throw new IOException("empty sample in "+line);
			if(pop.trim().isEmpty())  throw new IOException("empty sample in "+line);
			String prevpo= this.sample2population.get(sample);
			if(prevpo!=null && !prevpo.equals(pop))
				throw new IOException("two pop declared for "+sample);
			sample2population.put(sample, pop);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args, outputFile);
		}

	public static void main(String[] args)
		{
		new VcfGroupByPopulation().instanceMainWithExit(args);
		}
}
