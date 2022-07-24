/*

Copyright (c) 2022 Pierre Lindenbaum

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

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;
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
 java -jar dist/gff3upstreamorfasta  -R GRCh38.fa Homo_sapiens.GRCh38.107.chr.gff3.gz > uorf.gff3
```

note to self: test ENSG00000141736 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003529

END_DOC

*/
@Program(name="gff3upstreamorf",
description="Takes a standard GTF and generate a GTF containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"gff","gff3","uorf","uorf"},
creationDate="20220724",
modificationDate="20220724",
generate_doc=false
)
public class Gff3UpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(Gff3UpstreamOrf.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"--strength"},description="only accept events that are greater or equal to this Kozak strength.")
	private KozakSequence.Strength user_kozak_strength = KozakSequence.Strength.nil;
	@Parameter(names={"--break-original-orf"},description="if ATG(uORF) is in frame with original ORF , do not calculate the peptide beyond the original ATG.")
	private boolean break_in_original_orf = false;

	private static int ID_GENERATOR=0;
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private final GeneticCode geneticCode = GeneticCode.getStandard();

	
	
	private double kozakStrengthToScore(final KozakSequence.Strength f) {
		switch(f)
			{
			case nil: return 0;
			case Weak: return 10;
			case Moderate:  return 100;
			case Strong:return 1000;
			default: throw new IllegalStateException();
			}
		}
	
	
	
		

	
	
	
	
	/** an ORF in an UpstreamORF, a subsection of an mRNA */
	private class OpenReadingFrame extends AbstractCharSequence implements Locatable
		{
		private final int uniq_id;
		private final GFF3Transcript transcript;
		private KozakSequence kozak = null;
		private int _in_rna_atg0 = 0;
		private String peptide = null;
		//boolean uorf_atg_in_frame = false;
		//private final int chromStart0;
		OpenReadingFrame(final GFF3Transcript transcript) {
			this.transcript = transcript;
			this.uniq_id = (++ID_GENERATOR);
			}
		
		String getId() {
			return transcript.getId()+"."+this.uniq_id;
			}
		String getName() {
			return transcript.getName();
			}
		
		@Override
		public String getContig() {
			return transcript.getContig();
			}
		
		@Override
		public int getStart() {
			return this.transcript.convertRna0ToGenomicIndex1(
				getATGPos0InRNA() + (isPositiveStrand()?0:length()-1)
				);
			}
		
		@Override
		public int getEnd() {
			return this.transcript.convertRna0ToGenomicIndex1(
				getATGPos0InRNA() + (isPositiveStrand()?length()-1:0)
				);
			}
		
		boolean isPositiveStrand() {
			return transcript.isPositiveStrand();
			}
		
		public int getATGPos0InRNA() {
			return this._in_rna_atg0;
		}
		
		public int getATGPos1InGenome() {
			return this.transcript.convertRna0ToGenomicIndex1(getATGPos0InRNA());
		}
		
		@Override
		public int hashCode() {
			return Integer.hashCode(getATGPos0InRNA())*31 + peptide.hashCode();
			}
		
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof OpenReadingFrame)) return false;
			final OpenReadingFrame other = OpenReadingFrame.class.cast(obj);
			if(length()!=other.length()) return false;
			if(!peptide.equals(other.peptide)) return false;
			return this.getATGPos1InGenome()== other.getATGPos1InGenome();
			}
		
		@Override
		public int length() {
			return this.peptide.length()*3;
			}
		
		@Override
		public char charAt(final int index) {
			if(index>=this.length()) return 'N';// if stop is beyond uORF, in whole mRNA
			return this.transcript.charAt(this.getATGPos0InRNA()+index);
			}
		public int getFrameAt(final int genomic1) {
			final int rna0 = this.transcript.convertGenomicIndex1ToRna0(genomic1);
			return (rna0-this.getATGPos0InRNA())%3;
			}
		public void write(final Gff3Feature gene,Gff3Writer w) throws IOException {
			final int genomic_orf_start1 = this.getStart();
			final int genomic_orf_end1 = this.getEnd();
			
			
			
			final Map<String,List<String>> tr_attributes = new HashMap<>();
			tr_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("transcript:"+getId()));
			tr_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(getName()));
			tr_attributes.put("kozak", Collections.singletonList(String.valueOf(this.kozak)));
			tr_attributes.put("kozak-strength", Collections.singletonList(String.valueOf(this.kozak.getStrength())));
			tr_attributes.put("peptide", Collections.singletonList(this.peptide));
			tr_attributes.put("peptide-length", Collections.singletonList(String.valueOf(this.length())));
			tr_attributes.put("transcript_id", Collections.singletonList(String.valueOf(getId())));
			tr_attributes.put("original-atg-pos", Collections.singletonList(String.valueOf(transcript.getATGPos1InGenome())));
			tr_attributes.put("atg-pos", Collections.singletonList(String.valueOf(getATGPos1InGenome())));
			tr_attributes.put("in-frame-with-transcript-atg", Collections.singletonList(String.valueOf(this.getATGPos0InRNA()%3==this.transcript.getATGPos0InRNA()%3)));
			if(gene!=null) 	tr_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList(gene.getID()));

			final Gff3FeatureImpl trFeat = new Gff3FeatureImpl(
				getContig(),
				Gff3UpstreamOrf.class.getSimpleName(),
				"transcript",
				this.transcript.getStart(),
				this.transcript.getEnd(),
				kozakStrengthToScore(this.kozak.getStrength()),
				this.transcript.transcript.getStrand(),
				0,
				tr_attributes
				);
			w.addFeature(trFeat);
			
			int exon_idx=0;
			for(Gff3Feature exon: transcript.exons) {
				final Map<String,List<String>> ex_attributes = new HashMap<>();
				ex_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList(trFeat.getID()));
				ex_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("exon:"+getId()+".exon"+exon_idx));
				ex_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(
					getName()+".exon."+exon_idx)	
					);
				exon_idx++;
				
				
				final Gff3FeatureImpl exFeat = new Gff3FeatureImpl(
						this.getContig(),
						exon.getSource(),
						exon.getType(),
						exon.getStart(),
						exon.getEnd(),
						exon.getScore(),
						exon.getStrand(),
						exon.getPhase(),
						ex_attributes
						);
				w.addFeature(exFeat);
				}
			
			int cds_idx=0;
			for(Gff3Feature exon: transcript.exons) {
				if(exon.getEnd() < genomic_orf_start1 ) continue;
				if(exon.getStart() > genomic_orf_end1) continue;
				int cds_start = Math.max(genomic_orf_start1,exon.getStart());
				int cds_end = Math.min(genomic_orf_end1,exon.getEnd());
				
				final Map<String, List<String>> cds_attributes = new HashMap<>();
				cds_attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList(trFeat.getID()));
				cds_attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("cds:"+ transcript.getId()+".cds"+cds_idx));
				cds_attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(
					getName()+".cds."+(cds_idx+1))	
					);
				cds_idx++;

				
				
				final Gff3FeatureImpl cdsfeat = new Gff3FeatureImpl(
					getContig(),
					Gff3UpstreamOrf.class.getSimpleName(),
					"CDS",
					cds_start,
					cds_end,
					kozakStrengthToScore(this.kozak.getStrength()),
					transcript.transcript.getStrand(),
					getFrameAt(cds_start),
					cds_attributes
					);
				w.addFeature(cdsfeat);
				}
			}
		}
	
	/** mutated version of UpstreamORF */
	
	private boolean acceptKozak(final KozakSequence k) {
		return k.getStrength().compareTo(this.user_kozak_strength)<=0;
		}
	
	/**
	 * Wrapper around a GFF3 transcript
	 *
	 */
	private class GFF3Transcript extends AbstractCharSequence implements Locatable {
		final String fixedContig;
		final Gff3Feature transcript;
		final List<Gff3Feature> exons = new ArrayList<>();
		final List<Interval> CDS = new ArrayList<>();
		final int transcript_length;
		GFF3Transcript(final Gff3Feature transcript,final String fixedContig) {
			this.fixedContig = fixedContig;
			this.transcript = transcript;
			for(Gff3Feature c:transcript.getChildren()) {
				if(c.getType().equals("exon")) {
					exons.add(c);
					}
				else if(c.getType().equals("CDS")) {
					CDS.add(new Interval(c));
					}
				}
			Collections.sort(exons,(A,B)->Integer.compare(A.getStart(),B.getStart()));
			Collections.sort(CDS,(A,B)->Integer.compare(A.getStart(),B.getStart()));
			this.transcript_length = exons.stream().mapToInt(E->E.getLengthOnReference()).sum();
			}
		
		String getId() {
			String s = this.transcript.getID();
			if(s.startsWith("transcript:")) s=s.substring(11);
			return s+".uorf";
			}
		
		String getName() {
			final String s= this.transcript.getName();
			return StringUtils.isBlank(s)?getId():s;
			}
		
		@Override
		public String getContig() {
			return this.fixedContig;
			}
		@Override
		public int getStart() {
			return this.transcript.getStart();
			}
		
		@Override
		public int getEnd() {
			return this.transcript.getEnd();
			}
		
		boolean isPositiveStrand() {
			return transcript.getStrand() == Strand.POSITIVE;
			}
		
		int getATGPos1InGenome() {
			return isPositiveStrand()?
				CDS.get(0).getStart():
				CDS.get(CDS.size()-1).getEnd();
			}
		int getATGPos0InRNA() {
			return convertGenomicIndex1ToRna0(getATGPos1InGenome());
			}
		
		/** convert 0-based rna position to 1-based genomic position */
		public int convertRna0ToGenomicIndex1(int p0) {
			final int p0src = p0;
			if(isPositiveStrand()) {
				for(int i=0;i< this.exons.size();i++) {
					final Locatable exon = this.exons.get(i);
					final int exonLength = exon.getLengthOnReference();
					if(p0 >= exonLength ) {
						p0-= exonLength;
						continue;
						}
					return exon.getStart()+p0;
					}
				}
			else
				{
				for(int i=this.exons.size()-1;i>=0;i--) {
					final Locatable exon = this.exons.get(i);
					final int exonLength = exon.getLengthOnReference();
					if(p0 >= exonLength ) {
						p0-= exonLength;
						continue;
						}
					return exon.getEnd()-p0;
					}
				}
			throw new IndexOutOfBoundsException("p0="+p0src+" length:"+length());
			}
		/** convert 1-based genomic position to 0-based rna position */
		public int convertGenomicIndex1ToRna0(final int p1) {
			if(isPositiveStrand()) {
				int n=0;
				for(int i=0;i< this.exons.size();i++) {
					final Locatable exon = this.exons.get(i);
					if(exon.getStart() <= p1 && p1 <= exon.getEnd()) {
						return n+ (p1-exon.getStart()); 
						}
					n += exon.getLengthOnReference();
					}
				}
			else
				{
				int n=0;
				for(int i=this.exons.size()-1;i>=0;i--) {
					final Locatable exon = this.exons.get(i);
					if(exon.getStart() <= p1 && p1 <= exon.getEnd()) {
						return n+ (exon.getEnd()-p1); 
						}
					n += exon.getLengthOnReference();
					}
				}
			throw new IndexOutOfBoundsException("p0="+p1+" length:"+length());
			}
		@Override
		public char charAt(int index) {
			final char c= genomicSequence.charAt(convertRna0ToGenomicIndex1(index)-1);
			return Character.toUpperCase(isPositiveStrand()?c:AcidNucleics.complement(c));
			}
		
		@Override
		public int length() {
			return this.transcript_length;
			}
		final boolean containsOpenReadingFrame(final OpenReadingFrame orf) {
			final int atg1genome = orf.getATGPos1InGenome();
			if(this.getATGPos1InGenome()==atg1genome) return true;
			// other must have atg in frame and ATG before the observed one
			final int my_atg0 = getATGPos0InRNA();
			final int orf_atg0 = orf.getATGPos0InRNA();
			if(orf_atg0%3 != my_atg0%3) return false;
			if(orf_atg0 < my_atg0) return false;
			return atg1genome>= this.transcript.getStart() && atg1genome<=this.transcript.getEnd();
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
		
		public Set<OpenReadingFrame> getOpenReadingFrames() {
			final int transcript_atg0 = getATGPos0InRNA();
			final Set<OpenReadingFrame> orfs = new HashSet<>();

			for(int i=0;i+2< this.length() && i < transcript_atg0;i++) {
				if(!this.isATG(i)) continue;
				final KozakSequence kozak = new KozakSequence(this, i);
				if(!acceptKozak(kozak))  continue;
				
				
				
				OpenReadingFrame orf = new OpenReadingFrame(this);
				orf.kozak = kozak;
				orf._in_rna_atg0 = i;
				final StringBuilder pep = new StringBuilder();
				for(int j=i;j+2< this.length() ;j+=3) {
					if(j >= transcript_atg0 && j%3==transcript_atg0%3 && break_in_original_orf) break;
					
					final char aa = geneticCode.translate(
						charAt(j+0),	
						charAt(j+1),
						charAt(j+2)
						);
					if(aa=='?')  {/* may happen in large regions with NNNN */
						orf=null;
						break;
						}
					
					if(geneticCode.isStop(aa)) {
						break;
						}
					pep.append(aa);
					}
				if(orf==null || pep.length()==0) continue;
				orf.peptide = pep.toString();
				orfs.add(orf);
				}
			return orfs;
			}
		
		}
	
	
	private GFF3Transcript withTranscript(Gff3Feature transcript,final String fixedContig) {
		if(transcript.getStrand()==Strand.NONE) return null;
		final GFF3Transcript tr = new GFF3Transcript(transcript,fixedContig);
		for(final Locatable cds:tr.CDS) {
			if(!tr.exons.stream().anyMatch(E->E.contains(cds))) {
				LOG.error("in "+tr.getId()+" cds "+cds+" in not covered by some exon.");
				return null;
				}
			}
		
		return (tr.CDS.isEmpty() || tr.exons.isEmpty()) ? null: tr;
		}

	
	private void withGene(Gff3Feature gene,final ContigNameConverter ctgConverter, Gff3Writer w) throws IOException {
		final String ctg = ctgConverter.apply(gene.getContig());
		if(StringUtils.isBlank(ctg)) return ;
		if(genomicSequence==null || !genomicSequence.getContig().equals(ctg)) {
			genomicSequence= new GenomicSequence(this.indexedFastaSequenceFile, ctg);
			}
		ID_GENERATOR=0;
		final List<GFF3Transcript> transcripts = new ArrayList<>();
		final Set<OpenReadingFrame> orfs = new HashSet<>();
		for(Gff3Feature c:gene.getChildren()) {
			final GFF3Transcript tr = withTranscript(c,ctg);
			if(tr==null) continue;
			transcripts.add(tr);
			}
		if(transcripts.isEmpty()) return;
		
		
		for(GFF3Transcript tr:transcripts) {
			for(OpenReadingFrame orf: tr.getOpenReadingFrames()) {
				if(transcripts.stream().
					filter(TR->TR!=orf.transcript).	
					anyMatch(TR->TR.containsOpenReadingFrame(orf))) continue;
				orfs.add(orf);
				}
			}
		if(orfs.isEmpty()) return;
		
		final Map<String,List<String>> g_atts = new HashMap<>(gene.getAttributes());
		if(!gene.getID().startsWith("gene:")) {
			g_atts.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("gene:"+gene.getID()));
			}
	
		
		final Gff3FeatureImpl gFeat = new Gff3FeatureImpl(
				ctg,
				gene.getSource(),
				gene.getType(),
				gene.getStart(),
				gene.getEnd(),
				gene.getScore(),
				gene.getStrand(),
				gene.getPhase(),
				g_atts
				);
		w.addFeature(gFeat);
		for(OpenReadingFrame orf:orfs) {
			orf.write(gFeat,w);
			}
		}
	
	
	private void recursive(Gff3Feature feat,final ContigNameConverter ctgConverter ,Gff3Writer w)  throws IOException {
		if(feat.getType().equals("gene")) {
			withGene(feat,ctgConverter,w);
			}
		else
			{
			for(Gff3Feature c:feat.getChildren()) {
				recursive(c,ctgConverter,w);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {		
		try {
			this.indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(refDict);
			final String input = oneFileOrNull(args);
			
			htsjdk.tribble.readers.LineIterator lr =
						(StringUtils.isBlank(input)?
						IOUtils.openStreamForLineIterator(stdin()):
						IOUtils.openURIForLineIterator(input)
						);
						
			
			try(OutputStream ps = super.openPathOrStdoutAsStream(this.outputFile)) {
				try(Gff3Writer w = new Gff3Writer(ps){
					protected String escapeString(String s) {
						return s;
					}
				}) {
					
					for(final KozakSequence.Strength f:KozakSequence.Strength.values()) {
						w.addComment("kozak."+f.name()+"="+kozakStrengthToScore(f));
						}
					final String gtfSource = getProgramName().toLowerCase();
	
					for(final SAMSequenceRecord ssr: refDict.getSequences()) {
						w.addComment("contig "+ssr.getSequenceName()+": length:"+ssr.getSequenceLength());
					}
					
					w.addComment(gtfSource+":"+JVarkitVersion.getInstance().toString());
					
					if(!StringUtils.isBlank(input)) {
						w.addComment("gtf:"+input);
					}
	
					final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.DEEP,S->!(
						S.equals("biotype") ||
						S.equals(Gff3Constants.ID_ATTRIBUTE_KEY) ||
						S.equals(Gff3Constants.PARENT_ATTRIBUTE_KEY) ||
						S.equals(Gff3Constants.NAME_ATTRIBUTE_KEY)
						));
					codec.readHeader(lr);
					while(!codec.isDone(lr)) {
						final Gff3Feature feat=codec.decode(lr);
						if(feat==null) continue;
						recursive(feat,ctgConverter,w);
						}
					codec.close(lr);
					}
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
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		
		}
	
	public static void main(final String[] args) {
		new Gff3UpstreamOrf().instanceMainWithExit(args);
	}
}