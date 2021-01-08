/*

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
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.bio.structure.UTR;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
/**
BEGIN_DOC

## inspiration

Part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

Wikipedia:

> An Upstream Open Reading Frame (uORF) is an open reading frame (ORF) within the 5' untranslated region (5'UTR) of an mRNA. uORFs can regulate eukaryotic gene expression.
> Translation of the uORF typically inhibits downstream expression of the primary ORF. In bacteria, uORFs are called leader peptides, and were originally discovered on the basis of their impact on the regulation of genes involved in the synthesis or transport of amino acids. 

## Examples

### Example 1

```

```

note to self: test ENSG00000141736 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003529

END_DOC

*/
@Program(name="gtfupstreamorf",
description="Takes a standard GTF and generate a GTF containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"gtf","uorf"},
creationDate="20190718",
modificationDate="20200804",
generate_doc=false
)
public class GtfUpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(GtfUpstreamOrf.class).make();
	private static final int NPOS=-1;
	private static final Paranoid PARANOID = Paranoid.createThrowingInstance();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"--strength"},description="only accept events that are greater or equal to this Kozak strength.")
	private KozakSequence.Strength user_kozak_strength = KozakSequence.Strength.nil;


	
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;
	private final GeneticCode geneticCode = GeneticCode.getStandard();

	
	
	private final static char NO_BASE='\0';
	
	
	private int kozakStrengthToScore(final KozakSequence.Strength f) {
		switch(f)
			{
			case nil: return 0;
			case Weak: return 10;
			case Moderate:  return 100;
			case Strong:return 1000;
			default: throw new IllegalStateException();
			}
		}
	
	
	
	/** Sequence for a Transcript */
	private class RNASequence extends AbstractCharSequence {
		private final int mrnaIndex0ToGenomic0[];
		private final char mrnaIndex0ToBase[];
		private final Map<Integer, Integer> genomic0ToRNAindex0;
		private final Transcript transcript;

		RNASequence(final Transcript transcript) {
			this.transcript = transcript;
			int mrna_index=0;
			final int transcript_length = transcript.getTranscriptLength();
			this.mrnaIndex0ToGenomic0 = new int[transcript_length];
			this.mrnaIndex0ToBase = new char[transcript_length];
			Arrays.fill(this.mrnaIndex0ToBase,NO_BASE);
			Arrays.fill(this.mrnaIndex0ToGenomic0,-1);
			this.genomic0ToRNAindex0 = new HashMap<>(transcript_length);
			
			
			if(transcript.isPositiveStrand())
				{
				for(final Exon ex:this.transcript.getExons())
					{
					for(int pos1=ex.getStart();pos1<=ex.getEnd();++pos1)
						{
						this.mrnaIndex0ToGenomic0[mrna_index] =  pos1-1;
						this.genomic0ToRNAindex0.put( pos1-1,mrna_index);
						++mrna_index;
						}
					}
				}
			else
				{
				for(int i=transcript.getExonCount()-1;i>=0;i--)
					{
					final Exon ex = transcript.getExon(i);
					for(int pos1=ex.getEnd();pos1>=ex.getStart();--pos1)
						{
						this.mrnaIndex0ToGenomic0[mrna_index] = pos1-1;
						this.genomic0ToRNAindex0.put( pos1-1,mrna_index);
						++mrna_index;
						}
					}
				}
			if(mrna_index!=this.mrnaIndex0ToGenomic0.length) {
				throw new IllegalArgumentException("Cannot fill all genomic positions for "+transcript.getId()+" stop at "+mrna_index+" expected "+this.mrnaIndex0ToGenomic0.length);
				}
			}
		
		@Override
		public int length() {
			return mrnaIndex0ToGenomic0.length;
			}
		
		
		/** return the associated Transcript */
		public Transcript getTranscript() {
			return this.transcript;
			}
		
		@Override
		public char charAt(final int index_in_rna0) {
			if(index_in_rna0< 0 || index_in_rna0>=this.length()) throw new IndexOutOfBoundsException(""+index_in_rna0+" "+this.length());
			if(this.mrnaIndex0ToBase[index_in_rna0]!=NO_BASE) {
				return this.mrnaIndex0ToBase[index_in_rna0];
				}
			final int g0  = this.mrnaIndex0ToGenomic0[index_in_rna0];
			char c = Character.toUpperCase(GtfUpstreamOrf.this.genomicSequence.charAt(g0));
			if(getTranscript().isNegativeStrand()) {
				c=AcidNucleics.complement(c);
				}
			this.mrnaIndex0ToBase[index_in_rna0] = c;
			if(!AcidNucleics.isIUPAC(c)) throw new IllegalArgumentException("not ATGC: "+c+" index_in_rna="+index_in_rna0+" g0="+g0+" +="+getTranscript().isPositiveStrand());
			return c;
			}
		
		/** return this transcript 0-based position of ATG in RNA */
		int getATG0InRNA() {
			if(!this.getTranscript().hasCodonStartDefined()) throw new IllegalStateException();

			final int index0_atg ;
			if(this.getTranscript().isPositiveStrand()) {
				final int genomic_atg0 = this.getTranscript().getCodonStart().get().getStart() - 1;
				index0_atg = this.genomic0ToRNAindex0.get(genomic_atg0);
				}
			else
				{
				final int genomic_atg0 = this.getTranscript().getCodonStart().get().getEnd() - 1;
				index0_atg = this.genomic0ToRNAindex0.get(genomic_atg0);
				}
			
			if(index0_atg < 0 || index0_atg >= this.length()) {
				throw new IllegalStateException("cannot get pos of ATG for "+getTranscript());
			}
			if(!isATG(index0_atg)) {
				LOG.warn("not atg ? at "+index0_atg +" "+getTranscript().getId()+" "+getTranscript().getStrand());
				}
			return index0_atg;
			}
		
		/** find best OenReading frame in the region */
		public Set<OpenReadingFrame> getUpstreamOpenReadingFrames() {
			if(!this.getTranscript().hasCodonStartDefined()) throw new IllegalStateException();
			if(!this.getTranscript().getCodonStart().isPresent()) {
				return Collections.emptySet();
				}
			
			final int atg0_in_mrna = getATG0InRNA();

			final Set<OpenReadingFrame> set = new HashSet<>();
			for(int i=0;i +2< atg0_in_mrna;i++) {
				if(!isATG(i)) continue;
				final KozakSequence kozak = new KozakSequence(this, i);
				if(!acceptKozak(kozak))  continue;
				
				final OpenReadingFrame orf = new OpenReadingFrame(this);
				final StringBuilder pep = new StringBuilder();
				orf.in_rna_atg0 = i;
				orf.kozak = kozak;
				orf.uorf_atg_in_frame = (i%3 == atg0_in_mrna%3);
				
				// if it's in frame, we don't scan beyond mrna-atg, otherwise it's the whole transcript
				final int stop_here = orf.uorf_atg_in_frame?atg0_in_mrna:this.length();
				
				for(int j=i;j+2< stop_here;j+=3) {
					final char aa = geneticCode.translate(
						charAt(j+0),	
						charAt(j+1),
						charAt(j+2)
						);
					if(aa=='?')  {
						// it happends e.g ENSG00000168385 overlaps a large region with 'NNNN'
						LOG.warn("bad amino acid in genomic region "+
						this.getTranscript().getContig()+":"+(1+this.mrnaIndex0ToGenomic0[j+0])+" "+
							this.getTranscript().getId()+" j="+j+" "+
							charAt(j+0)+charAt(j+1)+charAt(j+2)+" "+
							(this.mrnaIndex0ToGenomic0[j+0]+1)+"~"+
							(this.mrnaIndex0ToGenomic0[j+1]+1)+"~"+
							(this.mrnaIndex0ToGenomic0[j+2]+1)
							);
						break;
						}
					pep.append(aa);
					
					if(geneticCode.isStop(aa)) {
						orf.in_rna_stop0=j;
						break;
						}
					}
				orf.peptide = pep.toString();
				set.add(orf);
				}
			return set;
			}
		public List<Interval> getCodonBlocks(int r0a,int r0b,int r0c) {
			final int g0;
			final int g1;
			final int g2;
			if(isPositiveStrand())
				{
				g0 = this.mrnaIndex0ToGenomic0[r0a];
				g1 = this.mrnaIndex0ToGenomic0[r0b];
				g2 = this.mrnaIndex0ToGenomic0[r0c];
				}
			else
				{
				g0 = this.mrnaIndex0ToGenomic0[r0c];
				g1 = this.mrnaIndex0ToGenomic0[r0b];
				g2 = this.mrnaIndex0ToGenomic0[r0a];
				}
			if(g0+2==g2) {
				return Collections.singletonList(
						new Interval(getTranscript().getContig(),g0+1,g2+1)
						);
				}
			else if(g0+1==g1)
				{
				// System.out.println("#GOT splicing in codon ! :-) "+getTranscript().getContig()+" "+getTranscript().getId()+" "+g0+" "+g1+" "+g2);
				return Arrays.asList(
					new Interval(getTranscript().getContig(),g0+1,g1+1),
					new Interval(getTranscript().getContig(),g2+1,g2+1)
					);
				}
			else if(g1+1==g2)
				{
				// System.out.println("#GOT splicing in codon !! :-) "+getTranscript().getContig()+" "+getTranscript().getId()+" "+g0+" "+g1+" "+g2);
				return Arrays.asList(
					new Interval(getTranscript().getContig(),g0+1,g0+1),
					new Interval(getTranscript().getContig(),g1+1,g2+1)
					);
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		
		/** return true if the is an AT at this position */
		boolean isATG(final int index) {
			return index>=0 && 
					index+2 < this.length() &&
					charAt(index+0)=='A' &&
					charAt(index+1)=='T' && 
					charAt(index+2)=='G'
					;
			}
		
		
		/** get the strand of the gene */
		protected final boolean isPositiveStrand() {
			return this.getTranscript().isPositiveStrand();
			}
		}
	
	
	
	/** an ORF in an UpstreamORF, a subsectiobn of an mRNA */
	private class OpenReadingFrame extends AbstractCharSequence 
		{
		private final RNASequence mRNA;
		private KozakSequence kozak = null;
		private int in_rna_atg0 = NPOS;
		private int in_rna_stop0= NPOS;
		private String peptide = null;
		boolean uorf_atg_in_frame = false;
		//private final int chromStart0;
		OpenReadingFrame(final RNASequence mRNA) {
			this.mRNA = mRNA;
			}
		
		public RNASequence getTranscriptSequence() {
			return this.mRNA;
			}
		
		public final Transcript getTranscript() {
			return this.getTranscriptSequence().getTranscript();
			}
		
		@Override
		public int length() {
			return this.in_rna_stop0 - this.in_rna_atg0;
			}
		
		@Override
		public char charAt(final int index) {
			if(index>=this.length()) return 'N';// if stop is beyond uORF, in whole mRNA
			return this.getTranscriptSequence().charAt(this.in_rna_atg0+index);
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof OpenReadingFrame)) return false;
			final OpenReadingFrame other= OpenReadingFrame.class.cast(obj);
			if(!this.getTranscript().getContig().equals(other.getTranscript().getContig())) return false;
			final int x1=  this.mRNA.mrnaIndex0ToGenomic0[this.in_rna_atg0];
			final int x2=  other.mRNA.mrnaIndex0ToGenomic0[other.in_rna_atg0];
			if(x1!=x2) return false;
			return true;
			}
		
		@Override
		public int hashCode() {
			return peptide.hashCode();
			}
		public int getFrameAt(final int geomic1) {
			final int rna0 = this.mRNA.genomic0ToRNAindex0.get(geomic1-1);
			return (rna0-this.in_rna_atg0)%3;
			}
		}
	
	/** mutated version of UpstreamORF */
	
	private boolean acceptKozak(final KozakSequence k) {
		return k.getStrength().compareTo(this.user_kozak_strength)<=0;
		}
	
	private String keyvalue(final String key,Object value) {
		return key.replace('-', '_')+" \""+value+"\" ; ";
	}
	
	@Override
	public int doWork(final List<String> args) {
		GtfReader gtfReader = null;
		PrintWriter pw = null;
		try {
			this.indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);
			final ContigDictComparator ctgDictComparator = new ContigDictComparator(refDict);
			
			final String input = oneFileOrNull(args);
			gtfReader = input==null?new GtfReader(stdin()):new GtfReader(input);
			gtfReader.setContigNameConverter(this.refCtgNameConverter);
			
			final List<Gene> genes = gtfReader.getAllGenes().
					stream().
					filter(G->G.hasStrand()).
					sorted((A,B)->{
				int i = ctgDictComparator.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return Integer.compare(A.getEnd(), B.getEnd());
				}).collect(Collectors.toList());
			
			gtfReader.close();gtfReader=null;

			
			pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			for(final KozakSequence.Strength f:KozakSequence.Strength.values()) {
				pw.println("#kozak."+f.name()+"="+kozakStrengthToScore(f));
			}
			final String gtfSource = getProgramName().toLowerCase();

			
			
			for(final SAMSequenceRecord ssr: refDict.getSequences()) {
				pw.println("##contig "+ssr.getSequenceName()+": length:"+ssr.getSequenceLength());
			}
			
			pw.println("#"+gtfSource+":"+JVarkitVersion.getInstance().toString());
			
			if(!StringUtils.isBlank(input)) {
				pw.println("#gtf:"+input);
			}
			
			final ProgressFactory.Watcher<Gene> progress = ProgressFactory.newInstance().dictionary(refDict).logger(LOG).build();
			for(final Gene gene:genes)
				{
				progress.apply(gene);
				
				/* new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(gene.getContig())) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, gene.getContig());
					}

				
				final List<RNASequence> rnas = gene.getTranscripts().
						stream().
						filter(T->T.isCoding()).
						filter(T->T.hasStrand()).
						filter(T->T.hasCodonStartDefined()).
						filter(T->T.getTranscriptUTR5().isPresent()).
						map(T->new RNASequence(T)).
						collect(Collectors.toList());

				if(rnas.isEmpty()) continue;
				
				final Set<OpenReadingFrame> orfs = rnas.
						stream().
						flatMap(R->R.getUpstreamOpenReadingFrames().stream()).
						collect(Collectors.toSet());

				if(orfs.isEmpty()) continue;
				
				boolean gene_printed = false;
				
				
				for(final OpenReadingFrame uORF : orfs)
					{
					
					/* is there any other RNA containing this uORF ?*/
					if(rnas.stream().
						//skip self
						filter(other->!other.getTranscript().getId().equals(uORF.mRNA.getTranscript().getId())).
						anyMatch(other->{
						// other must have atg in frame and ATG before the observed one
						final int other_atg0 = other.getATG0InRNA();
						if(uORF.in_rna_atg0%3 != other_atg0%3) return false;
						return uORF.in_rna_atg0 >= other_atg0;
						})) continue;
					
					final String transcript_id = uORF.getTranscript().getId()+".uorf"+(1+uORF.in_rna_atg0);
					
					if(!gene_printed) {
						gene_printed = true;
						pw.print(gene.getContig());
						pw.print("\t");
						pw.print(gtfSource);//source
						pw.print("\t");
						pw.print("gene");
						pw.print("\t");
						pw.print(gene.getStart());//start
						pw.print("\t");
						pw.print(gene.getEnd());//end
						pw.print("\t");
						pw.print(".");//score
						pw.print("\t");
						pw.print(gene.getStrand());//strand
						pw.print("\t");
						pw.print(".");//phase
						pw.print("\t");
						pw.print(keyvalue("gene_id",gene.getId()));
						pw.println();
						}

					
					
					// TRANSCRIPT
					final UTR utr_5_prime = uORF.getTranscript().getTranscriptUTR5().get();
					pw.print(gene.getContig());
					pw.print("\t");
					pw.print(gtfSource);
					pw.print("\t");
					pw.print("transcript");
					pw.print("\t");
					if(gene.isPositiveStrand()) {
						pw.print(uORF.mRNA.mrnaIndex0ToGenomic0[uORF.in_rna_atg0]+1);
						pw.print("\t");
						//pw.print(transcriptSequence.mrnaIndex0ToBase[uORF.in_rna_stop0]+1);
						if(uORF.in_rna_stop0==NPOS) {
							pw.print(uORF.mRNA.getTranscript().getTxEnd());
							}
						else
							{
							pw.print(uORF.mRNA.mrnaIndex0ToGenomic0[uORF.in_rna_stop0]+1);
							}
						} 
					else {
						if(uORF.in_rna_stop0==NPOS) {
							pw.print(uORF.mRNA.getTranscript().getTxStart());
							}
						else
							{
							pw.print(uORF.mRNA.mrnaIndex0ToGenomic0[uORF.in_rna_stop0]+1);
							}
						pw.print("\t");
						pw.print(uORF.mRNA.mrnaIndex0ToGenomic0[uORF.in_rna_atg0]+1);
						}
					pw.print("\t");
					pw.print(kozakStrengthToScore(uORF.kozak.getStrength()));//score
					pw.print("\t");
					pw.print(gene.getStrand());//strand
					pw.print("\t");
					pw.print("0");//phase
					pw.print("\t");
					pw.print(keyvalue("gene_id",gene.getId()));
					pw.print(keyvalue("transcript_id",transcript_id));
					pw.print(keyvalue("transcript_biotype","uORF"));
					pw.print(keyvalue("kozak-seq",uORF.kozak.getString()));
					pw.print(keyvalue("kozak-strength",uORF.kozak.getStrength()));
					pw.print(keyvalue("translation",uORF.peptide));
					pw.print(keyvalue("uORF-atg-in-frame-with-transcript-atg",uORF.uorf_atg_in_frame));
					pw.print(keyvalue("utr",utr_5_prime.toString()+" "+utr_5_prime.getStart()+"-"+utr_5_prime.getEnd()));
					
					pw.println();
					
					//Exon
					for(final Exon exon:uORF.getTranscript().getExons()) {
						pw.print(exon.getContig());
						pw.print("\t");
						pw.print(gtfSource);
						pw.print("\t");
						pw.print("exon");
						pw.print("\t");
						pw.print(exon.getStart());
						pw.print("\t");
						pw.print(exon.getEnd());
						pw.print("\t");
						pw.print(kozakStrengthToScore(uORF.kozak.getStrength()));//score
						pw.print("\t");
						pw.print(exon.getStrand());//strand
						pw.print("\t");
						pw.print(0);//phase
						pw.print("\t");
						pw.print(keyvalue("gene_id",gene.getId()));
						pw.print(keyvalue("transcript_id",transcript_id));
						pw.println();
					}
					
					
					
					final List<Interval> startBlocks = uORF.mRNA.getCodonBlocks(
							uORF.in_rna_atg0  ,
							uORF.in_rna_atg0+1,
							uORF.in_rna_atg0+2
							);
					
					final List<Interval> stopBlocks = uORF.in_rna_stop0!=NPOS?
								uORF.mRNA.getCodonBlocks(
								uORF.in_rna_stop0  ,
								uORF.in_rna_stop0+1,
								uORF.in_rna_stop0+2
								):Collections.emptyList();
					
					//CDS
					if(!stopBlocks.isEmpty()) {
						final int cdsStart = startBlocks.stream().mapToInt(B->B.getStart()).min().orElseThrow(IllegalStateException::new);
						final int cdsEnd = stopBlocks.stream().mapToInt(B->B.getEnd()).max().orElseThrow(IllegalStateException::new);
						for(final Exon exon:uORF.getTranscript().getExons()) {
							if(exon.getEnd() < cdsStart) continue;
							if(exon.getStart() > cdsEnd) break;
							pw.print(exon.getContig());
							pw.print("\t");
							pw.print(gtfSource);
							pw.print("\t");
							pw.print("CDS");
							pw.print("\t");
							pw.print(Math.max(cdsStart,exon.getStart()));
							pw.print("\t");
							pw.print(Math.min(cdsEnd,exon.getEnd()));
							pw.print("\t");
							pw.print(kozakStrengthToScore(uORF.kozak.getStrength()));//score
							pw.print("\t");
							pw.print(exon.getStrand());//strand
							pw.print("\t");
							pw.print(uORF.getFrameAt(Math.max(cdsStart,exon.getStart())));//phase
							pw.print("\t");
							pw.print(keyvalue("gene_id",gene.getId()));
							pw.print(keyvalue("transcript_id",transcript_id));
							pw.println();
						}
					}
								
					
					//CODON START
					for(final Interval startc: startBlocks) {
						PARANOID.assertLe(startc.getStart(), startc.getEnd());
						pw.print(startc.getContig());
						pw.print("\t");
						pw.print(gtfSource);
						pw.print("\t");
						pw.print("start_codon");
						pw.print("\t");
						pw.print(startc.getStart());
						pw.print("\t");
						pw.print(startc.getEnd());
						pw.print("\t");
						pw.print(kozakStrengthToScore(uORF.kozak.getStrength()));//score
						pw.print("\t");
						pw.print(gene.getStrand());//strand
						pw.print("\t");
						pw.print(uORF.getFrameAt(startc.getStart()));//phase
						pw.print("\t");
						pw.print(keyvalue("gene_id",gene.getId()));
						pw.print(keyvalue("transcript_id",transcript_id));
						pw.print(keyvalue("distance-mrna-atg",uORF.mRNA.getATG0InRNA()-uORF.in_rna_atg0));
						pw.print(keyvalue("pos0-in-mrna",uORF.in_rna_atg0));
						pw.print(keyvalue("spliced",startBlocks.size()>1));
						pw.println();
						}
					
					
					//CODON END
					for(final Interval stopc: stopBlocks) /* might be empty */ {
						PARANOID.assertLe(stopc.getStart(), stopc.getEnd());
						pw.print(stopc.getContig());
						pw.print("\t");
						pw.print(gtfSource);
						pw.print("\t");
						pw.print("stop_codon");
						pw.print("\t");
						pw.print(stopc.getStart());
						pw.print("\t");
						pw.print(stopc.getEnd());		
						pw.print("\t");
						pw.print(kozakStrengthToScore(uORF.kozak.getStrength()));//score
						pw.print("\t");
						pw.print(gene.getStrand());//strand
						pw.print("\t");
						pw.print(uORF.getFrameAt(stopc.getStart()));//phase
						pw.print("\t");
						pw.print(keyvalue("gene_id",gene.getId()));
						pw.print(keyvalue("transcript_id",transcript_id));
						pw.print(keyvalue("spliced",stopBlocks.size()>1));
						pw.println();
						}
					}
				
				}
			progress.close();
			
			pw.flush();
			pw.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		
		}
	
	public static void main(final String[] args) {
		new GtfUpstreamOrf().instanceMainWithExit(args);
	}
}