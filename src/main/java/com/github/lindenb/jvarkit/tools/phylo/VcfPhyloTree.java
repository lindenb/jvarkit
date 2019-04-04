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
package com.github.lindenb.jvarkit.tools.phylo;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
/**
BEGIN_DOC

## Example

```

```

END_DOC
 */
public class VcfPhyloTree extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfPhyloTree.class).make();
	private static int ID_GENERATOR=0;
	private static abstract class SampleNode
		{
		final int id= ++ID_GENERATOR;
		abstract void printDot(final PrintStream out);
		}
	
	private static class OneSampleNode extends SampleNode
		{
		private final String name;
		OneSampleNode(final String name)
			{
			this.name=name;
			}
		
		@Override void printDot(final PrintStream out)
			{
			out.println("n"+this.id+"[label=\""+this.name+"\"];");
			}
		@Override
		public String toString() {
			return name;
			}
		}
	private static class MergeSampleNode extends SampleNode
		{
		private final SampleNode n1;
		private final SampleNode n2;
		private final double distance;
		MergeSampleNode(final SampleNode n1,final SampleNode n2,double distance)
			{
			this.n1=n1;
			this.n2=n2;
			this.distance=distance;
			}
		
		@Override void printDot(final PrintStream out)
			{
			n1.printDot(out);
			n2.printDot(out);
			out.println("n"+this.id+"[shape=circle,label=\""+String.format("%.2f", this.distance)+"\"];");
			out.println("n"+this.n1.id+" -> n"+ this.id+";");
			out.println("n"+this.n2.id+" -> n"+ this.id+";");
			}

		
		@Override
		public String toString() {
			return "("+n1.toString()+","+n2+":"+this.distance+")";
			}
		}

	
	private static class GTPair
		{
		private final GenotypeType t1;
		private final GenotypeType t2;
		GTPair(final GenotypeType t1,final GenotypeType t2)
			{
			if(t1.compareTo(t2)<0) {
				this.t1 = t1;
				this.t2 = t2;
				}
			else
				{
				this.t1 = t2;
				this.t2 = t1;
				}
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + t1.hashCode();
			result = prime * result + t2.hashCode();
			return result;
			}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj)  return true;
			if (obj == null) {return false;}
			final GTPair other = (GTPair) obj;
			return this.t1.equals(other.t1) && this.t2.equals(other.t2);
			}
		
		}
	
	private final Map<GTPair,Double> gtype2score= new HashMap<>(); 
		
	private class Call
		{
		final List<GenotypeType> genotypes;
		Call(final Genotype g)
			{
			this.genotypes = Collections.singletonList(g.getType());
			}
		Call(final Call c1,final Call c2)
			{
			this.genotypes = new ArrayList<>(c1.genotypes.size()+c2.genotypes.size());
			this.genotypes.addAll(c1.genotypes);
			this.genotypes.addAll(c2.genotypes);
			}
		
		double distance(final Call c2) {
			int n=0;
			double d=0;
			for(final GenotypeType gt1:this.genotypes)
				{
				for(final GenotypeType gt2:c2.genotypes)
					{
					final GTPair gtpair = new GTPair(gt1, gt2);
					final Double score = VcfPhyloTree.this.gtype2score.get(gtpair);
					if(score ==null) throw new RuntimeException("diff score between "+gt1+" and "+gt2+" is not defined");
					d+=score.doubleValue();
					++n;
					}
				}
			return d/n;
			}
		}
	
	
	private VcfPhyloTree()
		{
		for(GenotypeType gt1:GenotypeType.values())
			{
			for(GenotypeType gt2:GenotypeType.values())
				{
				double score=5;
				if(gt1.equals(gt2))
					{
					score=0;
					}
				final  GTPair p=new GTPair(gt1, gt2);
				this.gtype2score.put(p, score);
				}
			}
		}
	

	@Override
	public int doWork(final List<String> args) {

		VCFIterator iter=null;
		try
			{
			iter = super.openVCFIterator(oneFileOrNull(args));
			
			final List<SampleNode> samplesNodes = new ArrayList<>(iter.getHeader().getSampleNamesInOrder().
					stream().map(S->new OneSampleNode(S)).
					collect(Collectors.toList()));
			

			if(samplesNodes.size()<3)
				{	
				LOG.error("expected at least 3 samples in input");
				return -1;
				}
			final List<List<Call>> variants = new ArrayList<>();
			while(iter.hasNext())
				{
				final VariantContext ctx=iter.next();
				final int n_nocall = (int) ctx.getGenotypes().stream().filter(G->G.isNoCall()).count();
				final int n_homref = (int) ctx.getGenotypes().stream().filter(G->G.isHomRef()).count();
				if(n_nocall + n_homref == ctx.getNSamples()) continue;
				
				
				final List<Call> row = new ArrayList<>(ctx.getNSamples());
				for(int i=0;i< samplesNodes.size();i++)
					{
					row.add(new Call(ctx.getGenotype(i)));
					}
				variants.add(row);
				}
			iter.close();
			iter=null;
			
			while(samplesNodes.size()>1)
				{
				LOG.info("Iteration N="+samplesNodes.size()+" / variants:"+variants.size());
				int best_x1 = -1;
				int best_x2 = -1;
				double best_distance=-1;
				for(int x1=0;x1+1 < samplesNodes.size();++x1)
					{
					for(int x2=x1+1;x2 < samplesNodes.size();++x2)
						{	
						double distance=0;
						for(final List<Call> row:  variants)
							{
							final Call c1 = row.get(x1);
							final Call c2 = row.get(x2);
							distance += c1.distance(c2);
							}
						if(best_distance==-1 || best_distance > distance)
							{
							best_distance=distance;
							best_x1 = x1;
							best_x2 = x2;
							}
						}
					}
				
				final List<SampleNode> newSampleNodes = new ArrayList<>(samplesNodes.size()-1);
				for(int x1=0;x1 < samplesNodes.size();++x1)
					{
					if(x1==best_x1 || x1==best_x2) continue;
					newSampleNodes.add(samplesNodes.get(x1));
					}
				
				LOG.info("Merging\n\t"+samplesNodes.get(best_x1)+"\nand\n\t"+samplesNodes.get(best_x2)+"\ndistance\n\t"+best_distance);

				
				newSampleNodes.add(new MergeSampleNode(
						samplesNodes.get(best_x1),
						samplesNodes.get(best_x2),
						best_distance)
						);
				
				for(int y=0;y< variants.size();++y)
					{
					final List<Call> row = variants.get(y);
					final List<Call> newrow = new ArrayList<>(samplesNodes.size()-1);
					for(int x1=0;x1 < samplesNodes.size();++x1)
						{
						if(x1==best_x1 || x1==best_x2) continue;
						newrow.add(row.get(x1));
						}
					newrow.add(new Call(row.get(best_x1),row.get(best_x2)));
					variants.set(y, newrow);
					}
				samplesNodes.clear();
				samplesNodes.addAll(newSampleNodes);
				
				}
			
			System.out.println("digraph G {");
			samplesNodes.get(0).printDot(System.out);
			System.out.println("}");
			
			
			LOG.info("Done");
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	public static void main(String[] args) {
		new VcfPhyloTree().instanceMainWithExit(args);
		}
	}
