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

*/
package com.github.lindenb.jvarkit.tools.skat;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherH;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenMAF;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

END_DOC

 */
@Program(
		name="vcfoptimizeped4skat",
		description="Optimize ped file for SKAT",
		keywords={"vcf","pedigree","skat","burden"}
		)
public class VcfOptimizePedForSkat extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfOptimizePedForSkat.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from  the VCFheader.")
	private File pedigreeFile=null;
	@Parameter(names={"-n","--remove"},description="max number of samples to remove")
	private int nSamplesRemove=1;
	@Parameter(names={"-seed","--seed"},description="random seed; -1=currentTimeMillis")
	private long seed=-0L;
	@Parameter(names={"--max-results"},description="max number of results.")
	private int max_results_output=10;
	@Parameter(names={"--max-iter"},description="max number of iterations. -1 == infinite")
	private long max_iterations=-1L;
	@Parameter(names={"--bootstrap"},description="bootstrap samples. Multiple list of sample separated with space, comma or semicolons")
	private String bootstrapSamples=null;
	
	@ParametersDelegate
	private SkatFactory skatInstance = new SkatFactory();

	
	private Pedigree pedigree=null;
	private Random random = null;
	private final List<Solution> bestSolutions = new ArrayList<>();
	private Double firstScore=null;
	
	/* collector for Variants after BurdenMaf & BurdenFIsherH */
	private static class VariantCollector implements VariantContextWriter
		{
		@SuppressWarnings("unused")
		VCFHeader header;
		final List<VariantContext> variants = new ArrayList<>();

		@Override
		public void setHeader(final VCFHeader header) {
			this.header = header;
			}

		@Override
		public void writeHeader(final VCFHeader header) {
			this.header = header;
			}
		@Override
		public void add(VariantContext vc) {
			if(vc.isFiltered()) return;
			this.variants.add(vc);
			}
		@Override public boolean checkError() {return false;}
		@Override public void close() {}
		}
	
	
	private class Solution implements Comparable<Solution>
		{
		String origin="";
		long generation;
		SkatFactory.SkatResult result;
		final Set<String> sampleSet= new TreeSet<>();
		
		
		@Override
		public int compareTo(final Solution o)
			{
			return this.result.compareTo(o.result);
			}
		@Override
		public int hashCode() {
			return sampleSet.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			return sampleSet.equals(Solution.class.cast(obj).sampleSet);
			}
		
		@Override
		public String toString()
			{
			return String.format("%6.3e",result.getPValue())+"\t"+
					(firstScore==null?"N/A":String.format("%6.1e",firstScore/result.getPValue()) )+
					"\t("+generation+")\t"+
					origin+"\t"+
					String.join(";",this.sampleSet);
			}
		}
	
	private void exec(
			final long generation,
			final VCFHeader header,
			final List<VariantContext> variants, 
			final List<Pedigree.Person> samples,
			final SkatFactory.SkatExecutor skatExecutor
			)
		{
		final Solution solution = new Solution();
		solution.generation = generation;
		String origin = "random";
		if(generation!=0)
			{
			int nRemove = 1+this.random.nextInt(this.nSamplesRemove);
			if(generation==1 && this.bootstrapSamples!=null)
				{
				origin="bootstrap";
				for(final String sample:this.bootstrapSamples.split("[; ,]"))
					{
					if(StringUtil.isBlank(sample)) continue;
					if(!samples.stream().anyMatch(S->S.getId().equals(sample))) {
						throw new JvarkitException.UserError("Sample "+sample+" not found in effective pedigree.");
						}
					LOG.info("bootstraping with "+sample);
					solution.sampleSet.add(sample);
					}
				}
			else if(generation%5==0 && !this.bestSolutions.isEmpty()) {
				int sol_index = this.random.nextInt(Math.min(this.max_results_output, this.bestSolutions.size()));
				final List<String> list =  new ArrayList<>(this.bestSolutions.get(sol_index).sampleSet);
				if(list.size()>1 && this.random.nextBoolean())
					{
					origin="best-minus-random";
					list.remove(this.random.nextInt(list.size()));
					}
				else if(list.size()<nRemove)
					{
					origin="best-plus-random";
					list.add( samples.get(this.random.nextInt(samples.size())).getId());
					}
				solution.sampleSet.addAll(list);
				}
			else if(generation%7==0 && this.bestSolutions.size()>2)
				{
				final Set<String> set=new HashSet<>(this.bestSolutions.get(0).sampleSet);
				set.addAll(this.bestSolutions.get(1).sampleSet);
				final List<String> bestsamples =  new ArrayList<>(set);
				Collections.shuffle(bestsamples, this.random);
				while(bestsamples.size()>nRemove) {
					bestsamples.remove(0);
					}
				solution.sampleSet.addAll(bestsamples);
				origin="best0-plus-best1";
				}
			else
				{
				while(nRemove>0)
					{
					final String sampleId ;
					
					if(generation%3==0L && 
							nRemove%2==0 && 
							this.bestSolutions.size()>0 && 
							!this.bestSolutions.get(0).sampleSet.isEmpty())
						{
						final List<String> bestsamples =  new ArrayList<>(this.bestSolutions.get(0).sampleSet);
						sampleId = bestsamples.get(this.random.nextInt(bestsamples.size()));
						origin="random-plus-best0";
						}
					else
						{
						sampleId = samples.get(this.random.nextInt(samples.size())).getId();
						}
					solution.sampleSet.add(sampleId);
					nRemove--;
					}
				}
			}
		else
			{
			origin="original";
			}
		if(generation>0 && solution.sampleSet.isEmpty()) return;
		while(solution.sampleSet.size()> this.nSamplesRemove)
			{
			LOG.warn("Hum... to many for "+origin);
			final List<String> L =  new ArrayList<>(solution.sampleSet);
			while(L.size()>0 && L.size()> this.nSamplesRemove) L.remove(this.random.nextInt(L.size()));
			solution.sampleSet.clear();
			solution.sampleSet.addAll(L);
			}
		if(this.bestSolutions.contains(solution)) return;
		solution.origin = origin;
		
		final List<Pedigree.Person> ped2 = new ArrayList<>(samples);
		ped2.removeIf(I->solution.sampleSet.contains(I.getId()));
		if(ped2.isEmpty()) return;
		if(!ped2.stream().anyMatch(P->P.isAffected())) return;
		if(!ped2.stream().anyMatch(P->P.isUnaffected())) return;
		
		// test MAF et Fisher
		
		final VCFHeader header2 = new VCFHeader(
				header.getMetaDataInInputOrder(),
				header.getSampleNamesInOrder()
				);
		final VcfBurdenFisherH.CtxWriterFactory fisherhFactory = new VcfBurdenFisherH.CtxWriterFactory();
		fisherhFactory.setCaseControlExtractor((H)->new HashSet<>(ped2));
		fisherhFactory.setIgnoreFiltered(true);
		final VcfBurdenMAF.CtxWriterFactory mafFactory = new VcfBurdenMAF.CtxWriterFactory();
		mafFactory.setCaseControlExtractor((H)->new HashSet<>(ped2));
		mafFactory.setIgnoreFiltered(true);

		
		final VariantCollector collector = new VariantCollector();
		VariantContextWriter vcw = mafFactory.open(collector);
		vcw = fisherhFactory.open(vcw);
		
		
		vcw.writeHeader(header2);
		for(final VariantContext vc: variants)
			{
			vcw.add(vc);
			}
		vcw.close();
		CloserUtil.close(fisherhFactory);
		CloserUtil.close(mafFactory);
		
		//
		solution.result = skatExecutor.execute(collector.variants, ped2);
		if(solution.result.isError()) return;
		if( this.bestSolutions.isEmpty() ||
			solution.compareTo(this.bestSolutions.get(this.bestSolutions.size()-1))<0)
			{
			this.bestSolutions.add(solution);
			if(this.firstScore==null)
				{
				this.firstScore = solution.result.getPValue();
				}
			
			Collections.sort(this.bestSolutions);
			final int BUFFER_RESULT=1000;
			while(this.bestSolutions.size()>BUFFER_RESULT) {
				this.bestSolutions.remove(this.bestSolutions.size()-1);
				}
			if(this.bestSolutions.indexOf(solution)<this.max_results_output) {
				
				final StringBuilder sb=new StringBuilder();
				sb.append(">>> ").append(generation).append("\n");
				
				this.bestSolutions.stream().
					limit(this.max_results_output).
					forEach(S->sb.append(S).append("\n"));
				sb.append("<<< ").append(generation).append("\n");
				
				if(outputFile==null)
					{
					stdout().println(sb.toString());
					}
				else
					{
					stderr().println(sb.toString());
					
					
					IOUtils.cat(sb.toString(),this.outputFile,false);
					}
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		VCFIterator r=null;
		try {
			this.random = new Random(this.seed==-1?System.currentTimeMillis():seed);
			if(this.nSamplesRemove<1) {
				LOG.error("bad number of samples to remove.");
				return -1;
			}
			
			final SkatFactory.SkatExecutor executor = this.skatInstance.build();
			
			r= super.openVCFIterator(oneFileOrNull(args));
			final List<VariantContext> variants = new ArrayList<>();
			while(r.hasNext())
				{
				final VariantContext ctx = r.next();
				if(!executor.getUpstreamVariantFilter().test(ctx)) continue;
				variants.add(ctx);
				}
			LOG.info("number of variants : "+variants.size());
			if(variants.stream().map(V->V.getContig()).collect(Collectors.toSet()).size()!=1) {
				LOG.error("multiple contig/chromosome in input.");
				return -1;
				}
			final VCFHeader header= r.getHeader();
			final Set<String> sampleNames = new HashSet<>(header.getSampleNamesInOrder());
			if(this.pedigreeFile!=null)
				{
				this.pedigree = new Pedigree.Parser().parse(this.pedigreeFile);
				}
			else
				{
				this.pedigree = new Pedigree.Parser().parse(header);
				}
			r.close();r=null;
			
			final List<Pedigree.Person> samples = this.pedigree.getPersons().
					stream().
					filter(P->(P.isAffected() || P.isUnaffected())).
					filter(S->sampleNames.contains(S.getId())).
					collect(Collectors.toList())
					;
			if(samples.isEmpty())
				{
				LOG.error("Not enough Samples in pedigree/vcf");
				return -1;
				}
			if(!samples.stream().anyMatch(P->P.isAffected())) 
				{
				LOG.error("No affected Samples in pedigree/vcf");
				return -1;
				}
			if(!samples.stream().anyMatch(P->P.isUnaffected())) 
				{
				LOG.error("No unaffected Samples in pedigree/vcf");
				return -1;
				}

			
			
			long nIter=0L;
			while(max_iterations==-1L || nIter<max_iterations)
				{
				exec(nIter,header,variants,samples,executor);
				++nIter;
				}
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	public static void main(String[] args)
		{
		new VcfOptimizePedForSkat().instanceMainWithExit(args);
		}
	}
