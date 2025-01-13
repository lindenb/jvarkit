/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.fishercombinatorics;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

END_DOC

 */
@Program(
		name="fishercnvgenealgo",
		description="Optimize gene pool for fisher-test and CNV vcf",
		keywords={"vcf","pedigree","cnv","fisher"},
		creationDate="20211102",
		modificationDate="20211102",
		generate_doc=false
		)
public class VcfFisherCombinatorics extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfFisherCombinatorics.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-ped","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile = null;
	@Parameter(names={"-bed","--bed"},description=BedLineReader.OPT_DESC+". 4th column contains the name of the interval.",required=true)
	private Path bedPath = null;
	@Parameter(names={"--max-af"},description="minimal allele frequency." + FractionConverter.OPT_DESC,converter=FractionConverter.class,hidden=false )
	private double max_af=1.0;
	@Parameter(names={"--max-genes"},description="Maximum number of gene in each solution")
	private int max_num_genes=100;
	@Parameter(names={"--bnd"},description="ignore variants with INFO/SVTYPE=BND")
	private boolean discard_bnd = false;
	@Parameter(names={"--filtered"},description="ignore FILTERed variants")
	private boolean discard_filtered = false;

	
	
	
	private static class GeneInfo implements Locatable,Comparable<GeneInfo> {
		String contig;
		String name;
		final List<Locatable> parts = new ArrayList<>();
		BitSet casesflags;
		BitSet controlsflags;
		public String getName() {
			return name;
			}
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return parts.get(0).getStart();
			}
		@Override
		public int getEnd() {
			return parts.get(parts.size()-1).getEnd();
			}
		
		@Override
		public int hashCode() {
			int i = getContig().hashCode();
			i = i*31  + Integer.hashCode(getStart());
			i = i*31  + Integer.hashCode(getEnd());
			i = i*31  + getName().hashCode();
			return i;
			}
		
		@Override
		public int compareTo(GeneInfo o) {
			int i= getContig().compareTo(o.getContig());
			if(i!=0) return i;
			i = Integer.compare(this.getStart(), o.getStart());
			if(i!=0) return i;
			i = Integer.compare(this.getEnd(), o.getEnd());
			if(i!=0) return i;
			i = getName().compareTo(o.getName());
			return i;
			}
		@Override
		public String toString() {
			return getName();
			}
		}
	
	private static class Solution implements Comparable<Solution>
		{
		long generation;
		final List<GeneInfo> genes = new ArrayList<>();
		Double result = null;
		
		@Override
		public int compareTo(final Solution o)
			{
			return this.result.compareTo(o.result);
			}
	
		
		@Override
		public String toString()
			{
			return "["+generation+"] "+result+"\t"+genes.stream().map(G->G.toString()).collect(Collectors.joining(" "));
			}
		}
	
	private Solution recursion(
		final List<GeneInfo> geneList,
		final List<Integer> casesIdx,
		final List<Integer> controlsIdx,
		final List<Integer> init_loop_indexes,
		final int max_level,
		Solution best
		) {
		if(init_loop_indexes.size()==max_level) {
			final Solution curr = new Solution();
			
			final BitSet cases_seen = new BitSet(casesIdx.size());
			final BitSet ctrls_seen = new BitSet(casesIdx.size());
			
			for(int x:init_loop_indexes) {
				final GeneInfo gi = geneList.get(x);
				cases_seen.or(gi.casesflags);
				ctrls_seen.or(gi.controlsflags);
			}
			
			
			int affected_alt = 0;
			int affected_ref = 0;
			for(int y=0;y< cases_seen.size();++y) {
				if(cases_seen.get(y)) {
					affected_alt++;
				} else
				{
					affected_ref = 0;
				}
			}
			int unaffected_alt = 0;
			int unaffected_ref = 0;
			for(int y=0;y< ctrls_seen.size();++y) {
				if(ctrls_seen.get(y)) {
					unaffected_alt++;
				} else
				{
					unaffected_ref = 0;
				}
			}

			curr.result = FisherExactTest.compute(
					affected_alt, affected_ref,
					unaffected_alt, unaffected_ref
					).getAsDouble();
			if(best==null || curr.compareTo(best)<0) {
				for(int x:init_loop_indexes) {
					curr.genes.add(geneList.get(x));
				}
				Collections.sort(curr.genes);
				System.out.println(curr);
				return curr;
				}
			else
				{
				return best;
				}
			}
		final List<Integer> loop_indexes = new ArrayList<>(init_loop_indexes);
		int offset = loop_indexes.isEmpty()?0:1 + loop_indexes.get(loop_indexes.size()-1);
		while(offset< geneList.size()) {
			boolean valid = true;
			for(int x=0;x< loop_indexes.size();++x) {
				if(geneList.get(x).overlaps(geneList.get(offset))) {
					valid = false;
					break;
					}
				}
			if(valid) {
				loop_indexes.add(offset);
				best = recursion(geneList,casesIdx,controlsIdx,loop_indexes,max_level,best);
				loop_indexes.remove(loop_indexes.size()-1);
				}
			++offset;
			}
		return best;
		}
	
	
	
	@Override
	public int doWork(final List<String> args) 	{
		try {
			
			final Pedigree pedigree = new PedigreeParser().parse(this.pedigreeFile);
			final IntervalTreeMap<List<GeneInfo>> geneInfoMap = new IntervalTreeMap<>();
			final Predicate<Genotype> isGtCnv = G->G.isHet() || G.isHomVar();

			
			try(VCFIterator vcfIter = super.openVCFIterator(super.oneFileOrNull(args))) {
				final VCFHeader header = vcfIter.getHeader();
				final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(header));
				if(!header.hasGenotypingData()) {
					LOG.error("No Genotype data in "+header);
					return -1;
					}
				final Map<String,Integer> sample2index = header.getSampleNameToOffset();

				final List<Integer> casesIdx = pedigree.getSamplesInVcfHeader(header).
						filter(S->S.isAffected()).
						map(S->S.getId()).
						map(S->sample2index.get(S)).
						sorted().
						collect(Collectors.toList());
				LOG.info("cases N="+casesIdx.size());
				if(casesIdx.isEmpty()) {
					LOG.error("No affected/cases sample in the input VCF");
					return -1;
					}

				final List<Integer> controlsIdx = pedigree.getSamplesInVcfHeader(header).
					filter(S->S.isUnaffected()).
					map(S->S.getId()).
					map(S->sample2index.get(S)).
					sorted().
					collect(Collectors.toList());
				LOG.info("controls N="+controlsIdx.size());
				if(controlsIdx.isEmpty()) {
					LOG.error("No unaffected/controls sample in the input VCF");
					return -1;
					}

				
				
				final Predicate<VariantContext> acceptCtx = V->{
					if(discard_bnd && V.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) return false;
					if(discard_filtered && V.isFiltered()) return false;
					if(max_af<1.0) {
						final Iterator<Integer> it = Stream.concat(controlsIdx.stream(),casesIdx.stream()).iterator();
						int ac=0;
						int an=0;
						while(it.hasNext()) {
							switch(V.getGenotype(it.next()).getType()) {
								case HOM_VAR: ac+=2; an+=2; break;
								case HOM_REF: ac+=0; an+=2; break;
								case HET: ac+=1; an+=2; break;
								default:break;
								}
							}
						if(an==0) return false;
						if(ac/(double)an > max_af) return false;
						}
					
					return true;
				};

				
				
				try(BedLineReader br = new BedLineReader(this.bedPath)) {
					
					while(br.hasNext()) {
						final BedLine bed = br.next();
						final String ctg = ctgConverter.apply(bed.getContig());
						if(StringUtil.isBlank(ctg)) continue;
						final GeneInfo geneInfo = new GeneInfo();
						geneInfo.casesflags = new BitSet(casesIdx.size());
						geneInfo.controlsflags = new BitSet(controlsIdx.size());
						geneInfo.contig = ctg;
						geneInfo.name = bed.getOrDefault(3, "undefined");
						geneInfo.parts.add(new SimpleInterval(bed).renameContig(ctg));
						
						final Interval key = new Interval(geneInfo);
						List<GeneInfo> L = geneInfoMap.get(key);
						if (L==null) {
							L = new ArrayList<>();
							geneInfoMap.put(key,L);
							}
						L.add(geneInfo);
						}
					}

				if(geneInfoMap.isEmpty()) {
					LOG.error("no gene found in "+this.bedPath);
					return -1;
					}
				LOG.info("reading variants...");
				while(vcfIter.hasNext()) {
					final VariantContext ctx = vcfIter.next();
					if(!acceptCtx.test(ctx)) continue;
					for(List<GeneInfo> gil: geneInfoMap.getOverlapping(ctx)) {
						for(GeneInfo gi:gil) {
							for(int y=0;y<casesIdx.size();++y) {
								final Genotype gt = ctx.getGenotype(casesIdx.get(y));
								if(!isGtCnv.test(gt)) continue;
								gi.casesflags.set(y);
								}
							for(int y=0;y<controlsIdx.size();++y) {
								final Genotype gt = ctx.getGenotype(controlsIdx.get(y));
								if(!isGtCnv.test(gt)) continue;
								gi.controlsflags.set(y);
								}
							}
						}
					} // end loop over variants
				
				/* remove Genes without variant , count sample per gene*/
				for(List<GeneInfo> gil: geneInfoMap.values()) {
					gil.removeIf(GI-> GI.casesflags.nextSetBit(0)==-1 && GI.controlsflags.nextSetBit(0)==-1);
				} // end remove
			
				/* load bed file */
				final List<GeneInfo> geneList = geneInfoMap.values().stream().
					filter(L->!L.isEmpty()).
					flatMap(L->L.stream()).
					sorted().
					collect(Collectors.toList());

				if(geneList.isEmpty()) {
					LOG.error("no gene found in "+this.bedPath);
					return -1;
					}
				LOG.info("N Genes="+geneList.size());
				
				Solution best=null;
				for(int i=2;i< Math.min(this.max_num_genes,geneList.size());i++) {
					LOG.info("sarting loop from "+i);
					best = recursion(geneList, casesIdx, controlsIdx, new ArrayList<>(), i, best);
					}
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					pw.println(best);
					pw.flush();
					}
				} // end vcfIn
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	public static void main(String[] args)
		{
		new VcfFisherCombinatorics().instanceMainWithExit(args);
		}
	}
