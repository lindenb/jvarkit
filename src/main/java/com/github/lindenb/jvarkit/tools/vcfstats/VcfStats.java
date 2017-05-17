/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/*

 */
@Program(name="vcfstats",
	description="VCF statitics",
	keywords={"vcf","stats"},
	terms=Term.ID_0000018
	)
public class VcfStats extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStats.class).make();

	@Parameter(names={"-o","--output"},description="output Directory or zip file",required=true)
	private File outputFile = null;
	
	@Parameter(names={"-kg","--knownGenes"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private File kgFile = null;
	
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	private Pedigree pedigree = null;

	
	@Parameter(names={"--prefix"},description="File/zip prefix")
	private String prefix = "tmp";
	
	@Parameter(names={"--select","-select"},description="Optional Jexl expression to use when selecting the adjacent variants")
    public ArrayList<String> selectExpressions = new ArrayList<String>();
    	
	private List<JexlVCMatchExp> jexls = new ArrayList<>();
	
	private ArchiveFactory archiveFactory=null;
	
	
	
	private static class PlotXY
		{
		void plot(double x,double y) {
			
			}	
		}

	
	private static class HistoGram
		{
		private final List<String> labels;
		private final Map<String,Integer> label2column= new LinkedHashMap<>();
		private final AutoMap<String,List<Integer>> rows;
		
		HistoGram(final List<String> labels)
			{
			this.labels = Collections.unmodifiableList(labels);
			for(int i=0;i< this.labels.size();++i) {
				this.label2column.put(this.labels.get(i),i);
				}
			this.rows=new AutoMap<>(new LinkedHashMap<>(),()->{
				final List<Integer> L = new ArrayList<>(labels.size());
				while(L.size()< labels.size()) L.add(0);
				return L;
				});
			}
		public List<Integer> get(final String lbl)
			{
			return this.rows.get(lbl);
			}
		public void incr (final String rowColumnY,final String columnLabelX)
			{
			Integer x= this.label2column.get(columnLabelX);
			if(x==null) throw new NoSuchElementException("Cannot find "+columnLabelX+" in "+this.labels);
			final List<Integer> row = get(rowColumnY);
			row.set(x, row.get(x)+1);
			}
		void print(final PrintWriter out)
			{
			out.println("Count"+this.labels.stream().collect(Collectors.joining("\t")));
			for(final String key:this.rows.keySet())
				{
				out.print(key);
				out.println(this.rows.get(key).stream().map(N->String.valueOf(N)).collect(Collectors.joining("\t")));
				}
			}
		}
	
	private class Stats
		{
		Counter<String> counter=new Counter<String>();
		Counter<Integer> alternate_alleles=new Counter<Integer>();
		Counter<Integer> depth=new Counter<Integer>();
		Counter<Integer> qual=new Counter<Integer>();
		Counter<SequenceOntologyTree.Term> snpEffSo=new Counter<SequenceOntologyTree.Term>();
		Counter<SequenceOntologyTree.Term> vepSo=new Counter<SequenceOntologyTree.Term>();

		
		private void watch(boolean is_in_coding,String prefix,Allele allele1,Allele allele2)
			{
			for(int i=0;i< 2;++i)
				{
				Character a1= simpleATGC(allele1);
				Character a2= simpleATGC(allele2);
				if(a1==null || a2==null)
					{
					if(allele1.getDisplayString().length()!=allele2.getDisplayString().length())
						{
						counter.incr(prefix+".indels");
						}
					}
				else
					{
					counter.incr(prefix+".substitutions");
					if(isTransition(a1, a2))
						{
						counter.incr(prefix+".transitions");
						}
					else if(isTransversion(a1, a2))
						{
						counter.incr(prefix+".transversions");
						}
					}
				if(!is_in_coding) break;
				prefix=prefix+".coding";
				}
			}
		
			
			
			
		
		}
	
	public VcfStats()
		{
		//this.selectExpressions.add("vc azd");
		}
	
	// https://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
	private static boolean isTransversion(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='C') return true;
		if(a1=='C' &&  a2=='A') return true;
		if(a1=='G' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='G') return true;
		return false;
		}

	private static boolean isTransition(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='G') return true;
		if(a1=='G' &&  a2=='A') return true;
		if(a1=='C' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='C') return true;
		return false;
		}

	
	
	private static Character simpleATGC(Allele al)
		{
		if(al==null) return null;
		String s=al.getBaseString().toUpperCase();
		if(s==null || s.equals(".") || s.length()!=1 ) return null;
		switch(s.charAt(0))
			{
			case 'A': case 'T': case 'G': case 'C': return s.charAt(0);
			default: return null;
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		
		try {
			if(this.pedigreeFile!=null)
				{
				this.pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			final Map<String,String> exprMap=new HashMap<>();
        	for(int i=0;i+1< this.selectExpressions.size();i+=2) {
        		exprMap.put(this.selectExpressions.get(i),this.selectExpressions.get(i+1));
        		}
        	this.jexls = VariantContextUtils.initializeMatchExps(exprMap);
				
			this.archiveFactory = ArchiveFactory.open(this.outputFile);
			scan(oneFileOrNull(args));
			archiveFactory.close();
			archiveFactory=null;
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		} finally
			{
			CloserUtil.close(archiveFactory);
			}
		
		}
	
	private void scan(String vcfInputStream) {
		VcfIterator iter = null;
		
		
		try
			{			
			iter= super.openVcfIterator(vcfInputStream);
			final VCFHeader header=iter.getHeader();
			
			final Set<String> intersectSamples = (
					this.pedigree==null?Collections.emptySet():
						this.pedigree.getPersons().stream().
						map(P->P.getId()).
						filter(S->header.getSampleNamesInOrder().contains(S)).
						collect(Collectors.toSet())
					);
			
			final Set<String> affectedSamples = 
					(this.pedigree==null?Collections.emptySet():
					this.pedigree.getPersons().stream().
						filter(P->P.isAffected()).
						map(P->P.getId()).
						filter(S->intersectSamples.contains(S)).
						collect(Collectors.toSet())
					);
			
			final Set<String> unaffectedSamples = 
					(this.pedigree==null?Collections.emptySet():
					this.pedigree.getPersons().stream().
						filter(P->P.isUnaffected()).
						map(P->P.getId()).
						filter(S->intersectSamples.contains(S)).
						collect(Collectors.toSet())
					);

			final Set<String> maleSamples = 
					(this.pedigree==null?Collections.emptySet():
					this.pedigree.getPersons().stream().
						filter(P->P.isMale()).
						map(P->P.getId()).
						filter(S->intersectSamples.contains(S)).
						collect(Collectors.toSet())
					);
			final Set<String> femaleSamples = 
					(this.pedigree==null?Collections.emptySet():
					this.pedigree.getPersons().stream().
						filter(P->P.isFemale()).
						map(P->P.getId()).
						filter(S->intersectSamples.contains(S)).
						collect(Collectors.toSet())
					);

			
			final HistoGram sample2type = new HistoGram(header.getSampleNamesInOrder());
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			while(iter.hasNext())
				{
				final VariantContext ctx=progress.watch(iter.next());
				if(header.hasGenotypingData())
					{
					
					/** iterator over variants */
					for(int i=0;i< ctx.getNSamples();++i) {
						final Genotype g = ctx.getGenotype(i);
						sample2type.incr(g.getType().name(), g.getSampleName());
						}
					
					
					final List<Genotype> altGenotypes = ctx.getGenotypes().
							stream().
							filter(G->!(G.isNoCall() || G.isHomRef() || G.isFiltered())).
							collect(Collectors.toList());
					
					}
				else
					
					{
					
					}
				
				}
			progress.finish();
			iter.close();
			iter=null;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			throw new RuntimeException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(vcfInputStream);
			}
		}		
			
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfStats().instanceMainWithExit(args);
		}
	}
