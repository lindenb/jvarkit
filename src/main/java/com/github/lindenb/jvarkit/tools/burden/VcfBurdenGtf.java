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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Intron;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.bio.structure.TranscriptInterval;
import com.github.lindenb.jvarkit.util.bio.structure.UTR;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFReader;

/**
 BEGIN_DOC
 
 
 END_DOC
 */

@Program(name="vcfburdengtf",
description="Run Burden On different part of transcripts",
keywords={"vcf","burden","gtf","case","control"},
creationDate="20190806",
modificationDate="20190809",
generate_doc=false
)
public class VcfBurdenGtf
extends AbstractVcfBurden
{
	private static final Logger LOG = Logger.build(VcfBurdenGtf.class).make();
	@Parameter(names={"-g","-gtf","--gtf"},description="GTF file",required=true)
	private Path gtfFile = null;
	@Parameter(names={"-f","--filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("vc.isSNP() && vc.getNAlleles()==2 && !vc.getFilters().contains(\"ZZ\")");
	@Parameter(names={"-t","--treshold"},description="fisher-test treshold. Discard results greater than this value.")
	private double fisherTreshold = 1e-5;
	@Parameter(names={"-intergenic","--intergenic"},description="Ignored if empty. If it is '*' run only on all intergenic regions found in the GTF. Otherwise, run only on intergenic regions on chromosome specified by the parameter.")
	private String intergenic_contig = null;

	
		
	private class SubPartOfTranscript implements Locatable {
		private List<Locatable> intervals = new ArrayList<>();
		private final Transcript transcript;
		private final String label;
		SubPartOfTranscript(final Transcript transcript,final String label) {
			this.transcript = transcript;
			this.label = label;
			}
		SubPartOfTranscript(final Transcript transcript,final String label,final List<Locatable> locs) {
			this(transcript,label);
			this.intervals.addAll(locs);
			}
		SubPartOfTranscript(final TranscriptInterval transcriptInterval) {
			this(transcriptInterval.getTranscript(),transcriptInterval.getName());
			this.intervals.add(transcriptInterval);
			}
		public Transcript getTranscript() {
			return this.transcript;
			}
		@Override
		public String getContig() {
			return getTranscript().getContig();
			}
		@Override
		public int getStart() {
			return this.intervals.get(0).getStart();
			}
		@Override
		public int getEnd() {
			return this.intervals.get(this.intervals.size()-1).getEnd();
			}
		}


	private boolean accept(final VariantContext ctx) {
		if(!variantFilter.test(ctx))return false;
		return true;
	}
	

	@Override
	protected void runBurden(PrintWriter pw, VCFReader vcfReader, VariantContextWriter vcw) throws IOException {
			final SAMSequenceDictionary vcfDict = SequenceDictionaryUtils.extractRequired(vcfReader.getHeader());
			final List<Gene> all_genes;
			try(GtfReader gtfReader = new GtfReader(this.gtfFile)){
				gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(vcfDict));
				all_genes =  gtfReader.getAllGenes().
						stream().
						filter(G->StringUtil.isBlank(this.intergenic_contig) || this.intergenic_contig.equals("*") || this.intergenic_contig.equals(G.getContig()) ).
						sorted(new ContigDictComparator(vcfDict).createLocatableComparator()).
						collect(Collectors.toCollection(ArrayList::new));
				}

			pw.print("#chrom");
			pw.print("\t");
			pw.print("start0");
			pw.print("\t");
			pw.print("end");
			pw.print("\t");
			pw.print("name");
			pw.print("\t");
			pw.print("length");
			pw.print("\t");
			pw.print("gene");
			pw.print("\t");
			pw.print("type");
			pw.print("\t");
			pw.print("strand");
			pw.print("\t");
			pw.print("transcript");
			pw.print("\t");
			pw.print("gene-id");
			pw.print("\t");
			pw.print("intervals");
			pw.print("\t");
			pw.print("p-value");
			pw.print("\t");
			pw.print("affected_alt");
			pw.print("\t");
			pw.print("affected_hom");
			pw.print("\t");
			pw.print("unaffected_alt");
			pw.print("\t");
			pw.print("unaffected_hom");
			pw.print("\t");
			pw.print("variants.count");
			pw.println();

			
			
			final List<SimpleInterval> all_intergenic = new ArrayList<>();
			
			if(!StringUtil.isBlank(this.intergenic_contig)) {
				for(final SAMSequenceRecord ssr: vcfDict.getSequences()) {
					if(!(this.intergenic_contig.equals("*") || this.intergenic_contig.equals(ssr.getSequenceName()))) continue;
					final BitSet filled = new BitSet(ssr.getSequenceLength()+2);
					
					all_genes.stream().
						filter(G->G.getContig().equals(ssr.getSequenceName())).
						forEach(G->filled.set(G.getStart(), 1 /* bit set is 0 based */ + Math.min(G.getEnd(),ssr.getSequenceLength())));
					int i=1;
					while(i< ssr.getSequenceLength()) {
						if(filled.get(i)) {
							i++;
							continue;
							}
						int j=i;
						while(j< ssr.getSequenceLength() && !filled.get(j)) {
							j++;
							}
						all_intergenic.add(new SimpleInterval(ssr.getSequenceName(),i,j));
						i=j+1;
						}
					all_genes.removeIf(G->G.getContig().equals(ssr.getSequenceName()));
					}
				all_genes.clear();
				}
			
			final ProgressFactory.Watcher<Gene> progress = ProgressFactory.newInstance().logger(LOG).dictionary(vcfDict).build();
			
			/* run genes */
			for(final Gene gene : all_genes) {
				progress.apply(gene);
				final IntervalTree<VariantContext> intervalTree = new IntervalTree<>();
				vcfReader.query(gene).stream().
					filter(V->accept(V)).
					forEach(V->intervalTree.put(V.getStart(), V.getEnd(), V));
				if(intervalTree.size()==0) continue;
				

				
				for(final Transcript transcript:gene.getTranscripts()) {
					final List<SubPartOfTranscript> parts = new ArrayList<>();
					parts.addAll(transcript.getExons().stream().map(R->new SubPartOfTranscript(R)).collect(Collectors.toList()));
					parts.addAll(transcript.getIntrons().stream().map(R->new SubPartOfTranscript(R)).collect(Collectors.toList()));
					
					final int intron_window_size=1000;
					final int intron_window_shift=500;
					for(final Intron intron:transcript.getIntrons()) {
						if(intron.getLengthOnReference()<=intron_window_size) continue;
						int start_pos=intron.getStart();
						while(start_pos + intron_window_size  <= intron.getEnd())
							{
							int xend = Math.min(intron.getEnd(), start_pos+intron_window_size-1);
							int xstart = xend-intron_window_size-1;
							parts.add(new SubPartOfTranscript(transcript,
									intron.getName()+".Sliding",
									Collections.singletonList(new SimpleInterval(intron.getContig(),xstart,xend)))
									);
							start_pos +=intron_window_shift;
							}
						}
					
					
					for(final UTR utr:transcript.getUTRs()) {
						parts.add(new SubPartOfTranscript(transcript,utr.getName(),utr.getIntervals()));
						}
					
					if(transcript.getExonCount()>1) {
						parts.add(new SubPartOfTranscript(
								transcript, 
								"AllExons",
								transcript.getExons().stream().map(E->E.toInterval()).collect(Collectors.toList())
								));
						}
					if(transcript.hasCodonStartDefined() &&
						transcript.hasCodonStopDefined() &&
						transcript.getAllCds().size()>1) {
						parts.add(new SubPartOfTranscript(
								transcript, 
								"AllCds",
								transcript.getAllCds().stream().map(E->E.toInterval()).collect(Collectors.toList())
								));
						}
					
					final int L = transcript.getTranscriptLength();
					final int index2genomic[] = new int[L];
					int pos=0;
					for(final Exon exon : transcript.getExons())
						{
						for(int i=exon.getStart();i<=exon.getEnd();i++) {
							index2genomic[pos]=i;
							pos++;
							}
						}
					final int window_size = 200;
					final int window_shift = 100;

					int array_index=0;
					while(array_index< index2genomic.length) {
						final List<Locatable> intervals = new ArrayList<>();
						int prev_pos=-1;
						int start_pos=index2genomic[array_index];
						int i=0;
						while(i< window_size && array_index+i < index2genomic.length)
							{
							final int curr_pos = index2genomic[array_index+i];
							if(i>0 && prev_pos+1!= curr_pos) {
								intervals.add(new SimpleInterval(transcript.getContig(),start_pos,prev_pos));
								start_pos = curr_pos;
								}
							prev_pos = curr_pos;
							i++;
							}
						intervals.add(new SimpleInterval(transcript.getContig(),start_pos,prev_pos));
						parts.add(new SubPartOfTranscript(transcript,"Sliding",intervals));
						array_index+=window_shift;
						}
					
					for(final SubPartOfTranscript part: parts) {
						final List<VariantContext> variants = new ArrayList<>();
						for(final Locatable loc: part.intervals) {
							Iterator<IntervalTree.Node<VariantContext>> iter= intervalTree.overlappers(loc.getStart(),loc.getEnd());
							while(iter.hasNext()) variants.add(iter.next().getValue());
							}
						if(variants.isEmpty()) continue;

						final FisherResult fisher = runFisher(variants);
						if(fisher.p_value> this.fisherTreshold) continue;
						
						if(vcw!=null) {
							for(final VariantContext ctx:variants) {
								vcw.add(new VariantContextBuilder(ctx).
										attribute(BURDEN_KEY, VCFUtils.escapeInfoField(part.label)).
										make());
							}
						}
						
						pw.print(part.getContig());
						pw.print("\t");
						pw.print(part.getStart()-1);
						pw.print("\t");
						pw.print(part.getEnd());
						pw.print("\t");
						pw.print(part.label);
						pw.print("\t");
						pw.print(part.getLengthOnReference());
						pw.print("\t");
						pw.print(transcript.getProperties().getOrDefault("gene_name", "."));
						pw.print("\t");
						pw.print(transcript.getProperties().getOrDefault("transcript_type", "."));
						pw.print("\t");
						pw.print(gene.getStrand());
						pw.print("\t");
						pw.print(transcript.getId());
						pw.print("\t");
						pw.print(gene.getId());
						pw.print("\t");
						pw.print(part.intervals.stream().map(R->String.valueOf(R.getStart())+"-"+R.getEnd()).collect(Collectors.joining(";")));
						pw.print("\t");
						pw.print(fisher.p_value);
						pw.print("\t");
						pw.print(fisher.affected_alt);
						pw.print("\t");
						pw.print(fisher.affected_hom);
						pw.print("\t");
						pw.print(fisher.unaffected_alt);
						pw.print("\t");
						pw.print(fisher.unaffected_hom);
						pw.print("\t");
						pw.print(variants.size());
						pw.println();
						}
					
					}
				}
			progress.close();
			
			final ProgressFactory.Watcher<SimpleInterval> progress2 = ProgressFactory.newInstance().logger(LOG).dictionary(vcfDict).build();

			/** scan intergenics ... */
			for(final SimpleInterval intergenic : all_intergenic) {
				progress2.apply(intergenic);
				final int intergenic_window_size=2000;
				final int intergenic_window_shifr=100;
				final List<SimpleInterval> parts = new ArrayList<>();

				if(intergenic.getLengthOnReference()<=intergenic_window_size) continue;
				int start_pos=intergenic.getStart();
				while(start_pos + intergenic_window_size  <= intergenic.getEnd())
					{
					int xend = Math.min(intergenic.getEnd(), start_pos+intergenic_window_size-1);
					int xstart = xend-intergenic_window_size-1;
					parts.add(new SimpleInterval(intergenic.getContig(),xstart,xend));
					start_pos +=intergenic_window_shifr;
					}
					
					
				for(final SimpleInterval part: parts) {
					final List<VariantContext> variants = vcfReader.query(part).
							stream().
							filter(V->accept(V)).
							collect(Collectors.toList());
					
					if(variants.isEmpty()) continue;
	
					final FisherResult fisher = runFisher(variants);
					if(fisher.p_value> this.fisherTreshold) continue;
					
					final String label = "intergenic_"+part.getStart()+"_"+part.getEnd();
					if(vcw!=null) {
						for(final VariantContext ctx:variants) {
							vcw.add(new VariantContextBuilder(ctx).
									attribute(BURDEN_KEY, VCFUtils.escapeInfoField(label)).
									make());
						}
					}
					
					pw.print(part.getContig());
					pw.print("\t");
					pw.print(part.getStart()-1);
					pw.print("\t");
					pw.print(part.getEnd());
					pw.print("\t");
					pw.print(label);
					pw.print("\t");
					pw.print(part.getLengthOnReference());
					pw.print("\t");
					pw.print(".");
					pw.print("\t");
					pw.print("intergenic");
					pw.print("\t");
					pw.print(".");
					pw.print("\t");
					pw.print(".");
					pw.print("\t");
					pw.print(".");
					pw.print("\t");
					pw.print(""+part.getStart()+"-"+part.getEnd());
					pw.print("\t");
					pw.print(fisher.p_value);
					pw.print("\t");
					pw.print(fisher.affected_alt);
					pw.print("\t");
					pw.print(fisher.affected_hom);
					pw.print("\t");
					pw.print(fisher.unaffected_alt);
					pw.print("\t");
					pw.print(fisher.unaffected_hom);
					pw.print("\t");
					pw.print(variants.size());
					pw.println();
					}
					
				}
			progress2.close();
	
		}
	
public static void main(final String[] args) {
	new VcfBurdenGtf().instanceMainWithExit(args);
	}
}
