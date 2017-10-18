/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
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
	
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from  the VCFheader.")
	private File pedigreeFile=null;
	@Parameter(names={"-n","--remove"},description="max number of samples to remove")
	private int nSamplesRemove=1;
	@Parameter(names={"-seed","--seed"},description="random seed; -1=currentTimeMillis")
	private long seed=-1l;
	@Parameter(names={"--max-results"},description="max number of results.")
	private int max_results=10;
	@Parameter(names={"--max-iter"},description="max number of iterations. -1 == infinite")
	private long max_iterations=-1L;

	
	private Pedigree pedigree=null;
	private Random random = null;
	private final List<Solution> bestSolutions = new ArrayList<>();
	
	private static class Solution implements Comparable<Solution>
		{
		Skat.SkatResult result;
		Set<String> samples= new TreeSet<>();
		@Override
		public int compareTo(final Solution o)
			{
			return this.result.compareTo(o.result);
			}
		@Override
		public String toString()
			{
			return String.valueOf(result.getPValue())+"\t"+String.join(";",this.samples);
			}
		}
	
	private void exec(
			final List<VariantContext> variants, 
			final List<Pedigree.Person> samples
			)
		{
		final Solution solution = new Solution();
		final List<Pedigree.Person> ped2 = new ArrayList<>(samples);
		int nRemove = 1+this.random.nextInt(this.nSamplesRemove);
		while(nRemove>0 && ped2.size()>2)
			{
			final String sampleId = ped2.remove(this.random.nextInt(ped2.size())).getId();
			solution.samples.add(sampleId);
			nRemove--;
			}
		if(!ped2.stream().anyMatch(P->P.isAffected())) return;
		if(!ped2.stream().anyMatch(P->P.isUnaffected())) return;
		final Skat skat = new Skat();
		solution.result = skat.execute(variants, ped2);
		if(solution.result.isError()) return;
		if( this.bestSolutions.isEmpty() ||
			solution.compareTo(this.bestSolutions.get(this.bestSolutions.size()-1))<0)
			{
			this.bestSolutions.add(solution);
			Collections.sort(this.bestSolutions);
			while(this.bestSolutions.size()>this.max_results) {
				this.bestSolutions.remove(this.bestSolutions.size()-1);
				}
			stdout().println(">>>");
			this.bestSolutions.stream().forEach(S->stdout().println(S));
			stdout().println("<<<\n");
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		VcfIterator r=null;
		try {
			this.random = new Random(this.seed==-1?System.currentTimeMillis():seed);
			if(this.nSamplesRemove<1) {
				LOG.error("bad number of samples to remove.");
				return -1;
			}
			
			r= super.openVcfIterator(oneFileOrNull(args));
			final List<VariantContext> variants = new ArrayList<>();
			while(r.hasNext())
				{
				final VariantContext ctx = r.next();
				if(ctx.isFiltered()) continue;
				if(ctx.getNAlleles()!=2) continue;
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
				exec(variants,samples);
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
