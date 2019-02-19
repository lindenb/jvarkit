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
package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## inspiration

part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

## Example

```
 wget -q -O - "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr1.vcf.bgz" |\
 bcftools annotate -x "INFO,FILTER" |\
 java -jar /home/lindenb/src/jvarkit-git/dist/vcfscanupstreamorf.jar \
 	-R human_g1k_v37.fasta  --uorf-only  --canonical  

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	89333	rs1008713359	A	G	283.15	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89334|cap-atg:2708|atg-cds:40|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGC|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89318|atg-stop:16|pep:.
1	89359	rs1327179626	C	T	3839.47	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:89359|cap-atg:2683|atg-cds:65|atg-frame:atg-out-of-cds-frame|kozak-seq:TGCATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89356|atg-stop:3|pep:M
1	89391	rs1332733110	T	C	2045.80	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89391|cap-atg:2651|atg-cds:97|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89382|atg-stop:9|pep:VNK
1	89555	rs1200434471	A	C	283.62	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89557|cap-atg:2485|atg-cds:263|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89452|atg-stop:105|pep:MKSQNVSQKIIYNVCVRKRQYPSNFESLHQKENSK,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89557|cap-atg:1312|atg-cds:7|atg-frame:atg-out-of-cds-frame|kozak-seq:GAAATGA|kozak-strength:Moderate|stop-frame:not-in-frame-stop|stop-pos:89556|atg-stop:1|pep:.
1	89560	rs1234719556	C	A	448.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:T|atg-pos:89562|cap-atg:2480|atg-cds:268|atg-frame:atg-out-of-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89541|atg-stop:21|pep:IKLKVKM,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:T|atg-pos:89562|cap-atg:1307|atg-cds:12|atg-frame:atg-in-cds-frame|kozak-seq:TGAATGA|kozak-strength:Weak|stop-frame:not-in-frame-stop|stop-pos:89555|atg-stop:7|pep:.
1	89624	rs1166058274	T	C	606.05	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89624|cap-atg:2418|atg-cds:330|atg-frame:atg-in-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89624|cap-atg:1245|atg-cds:74|atg-frame:atg-out-of-cds-frame|kozak-seq:ACAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89609|atg-stop:15|pep:VKELF
1	89718	rs865856422	A	G	1466.11	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89719|cap-atg:2323|atg-cds:425|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89719|cap-atg:1150|atg-cds:169|atg-frame:atg-out-of-cds-frame|kozak-seq:AAAATGA|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89710|atg-stop:9|pep:TKL
1	89831	rs1209426147	A	G	372.62	.	UORF_DEL_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:C|atg-pos:89832|cap-atg:2210|atg-cds:538|atg-frame:atg-out-of-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:C|atg-pos:89832|cap-atg:1037|atg-cds:282|atg-frame:atg-in-cds-frame|kozak-seq:CTTATGT|kozak-strength:Weak|stop-frame:in-frame-stop|stop-pos:89811|atg-stop:21|pep:TFAIYHT
1	89945	rs1376722481	G	C	297.51	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:G|atg-pos:89947|cap-atg:2095|atg-cds:653|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:G|atg-pos:89947|cap-atg:922|atg-cds:397|atg-frame:atg-out-of-cds-frame|kozak-seq:AATATGC|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89803|atg-stop:144|pep:MPLASVSHLAKPRLRSGKMEAISSWERRQRRWEYYVATYVCNLPYLAL
1	90032	rs866094671	C	T	14378.50	.	UORF_ADD_ATG=transcript:ENST00000466430.1|strand:-|utr-start:89295|utr-end:120932|alt:A|atg-pos:90032|cap-atg:2010|atg-cds:738|atg-frame:atg-in-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC,transcript:ENST00000495576.1|strand:-|utr-start:89551|utr-end:91105|alt:A|atg-pos:90032|cap-atg:837|atg-cds:482|atg-frame:atg-out-of-cds-frame|kozak-seq:TTCATGG|kozak-strength:Moderate|stop-frame:in-frame-stop|stop-pos:89951|atg-stop:81|pep:MGQLVSRAARETKPQCTFYSLCAHQTC
```


END_DOC

*/
@Program(name="vcfscanupstreamorf",
description="Scan BAM for upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs ",
keywords={"vcf","uorf"},
creationDate="2019-02-18"
)
public class VcfScanUpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfScanUpstreamOrf.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx = null;
	@Parameter(names={"-k","-K","--kg","-kg"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUri = KnownGene.getDefaultUri();
	@Parameter(names={"--canonical"},description="reduce the number of transcripts. Keep one if some share the same UTR")
	private boolean canonical_utr = false;
	@Parameter(names={"--uorf-only"},description="only print variants having something to say about an uorf")
	private boolean print_uorf_only = false;
	@Parameter(names={"--plus"},description="only transcripts on '+' strand. For debugging only",hidden=true)
	private boolean plus_strand_only = false;
	@Parameter(names={"--dac"},description="disable scan for ATG creation")
	private boolean disable_atg_create = false;
	@Parameter(names={"--dad"},description="disable scan for ATG deletion")
	private boolean disable_atg_delete = false;
	@Parameter(names={"--dsc"},description="disable scan for STOP creation")
	private boolean disable_stop_create = false;
	@Parameter(names={"--dsd"},description="disable scan for STOP deletion")
	private boolean disable_stop_delete = false;
	@Parameter(names={"--kal"},description="disable scan for Kozak change")
	private boolean disable_kozak_alteration= false;
	@Parameter(names={"--strong"},description="only accept events that are related to 'Strong' Kozak pattern.")
	private boolean only_strong_kozak= false;

	
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;
	private final IntervalTreeMap<List<KnownGene>> knownGeneMap = new IntervalTreeMap<>();
	private final GeneticCode geneticCode = GeneticCode.getStandard();
	

	private final VCFInfoHeaderLine infoAddATG = new VCFInfoHeaderLine(
			"UORF_ADD_ATG",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about de-novo ATG creation"
			);
	private final VCFInfoHeaderLine infoDelATG = new VCFInfoHeaderLine(
			"UORF_DEL_ATG",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about ATG deletion"
			);
	private final VCFInfoHeaderLine infoAddSTOP = new VCFInfoHeaderLine(
			"UORF_ADD_STOP",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about de-novo STOP creation"
			);
	private final VCFInfoHeaderLine infoDelSTOP = new VCFInfoHeaderLine(
			"UORF_DEL_STOP",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about STOP deletion"
			);
	
	private final VCFInfoHeaderLine infoKozakAlteration = new VCFInfoHeaderLine(
			"UORF_KOZAK_ALTERATION",
			VCFHeaderLineCount.UNBOUNDED,
			VCFHeaderLineType.String,
			"Information about kozak sequences alteration"
			);

	
	enum KozakStrength {
		Strong,Moderate,Weak,nil
	}
	
	/* kozak consensus length */
	private static final int KOZAK_LENGTH=7;
	/* position of ATG in the kozak sequence */
	private static final int KOZAK_ATG=3;
	
	/* a Kozak consensus */
	private class KozakSequence
		extends AbstractCharSequence
		{
		final CharSequence delegate;
		final int offset_beg;
		KozakSequence(final  CharSequence delegate,final int offset_atg) {
			this.delegate = delegate;
			this.offset_beg = offset_atg-KOZAK_ATG;
			}
		
		@Override
		public char charAt(int index) {
			final int idx= this.offset_beg + index;
			return idx >=0 && idx <delegate.length()?delegate.charAt(idx):'N';
			}
		
		private boolean hasATG() {
			if(charAt(KOZAK_ATG  )!='A') return false;
			if(charAt(KOZAK_ATG+1)!='T') return false;
			if(charAt(KOZAK_ATG+2)!='G') return false;
			return true;
			}

		
		KozakStrength getStrength() {
			char c0 = charAt(0);
			if(c0=='N') return KozakStrength.nil;
			char c6 = charAt(6);
			if(c6=='N') return KozakStrength.nil;
			if(!hasATG()) return KozakStrength.nil;
			
			
			boolean ok_chart_at_0 = c0 == 'A' || c0 == 'G' ;
			
  			if ( ok_chart_at_0 && c6 == 'G'){
  				return  KozakStrength.Strong;
  				}
  			else if (ok_chart_at_0 || c6 == 'G'){
  				return  KozakStrength.Moderate;
  				}
  			return KozakStrength.Weak;
			}
		
		@Override
		public int length() {
			return KOZAK_LENGTH;
			}
		
		void annotate(final StringBuilder sb) {
			sb.append("kozak-seq:").append(this.toString()).
			append("|kozak-strength:").
			append(this.getStrength().name()).
			append("|");
			}
		
		}
	
	private abstract class AbstractUTRSequence extends AbstractCharSequence
		implements Locatable
		{
		abstract KnownGene getGene();
		abstract protected int convertOrfToGenomicIndex0(int index);

		boolean isATG(final int index) {
			return index>=0 && 
					index+2 < this.length() &&
					charAt(index+0)=='A' &&
					charAt(index+1)=='T' && 
					charAt(index+2)=='G'
					;
			}
		/** return true if translation at index is stop */
		boolean isStop(final int index) {
			return index>=0 && 
					index+2 < this.length() && 
					VcfScanUpstreamOrf.this.geneticCode.isStopCodon(
							charAt(index),
							charAt(index+1),
							charAt(index+2)
							)
					;
			}
		
		/** find best ATG upstream of stop 'best' is first in frame with atg with best kozak. returns -1 if not found*/
		int findBestATG(final int stop_pos) {
			int best=-1;
			for(int idx=stop_pos-1;idx>=0;--idx) {
				if(idx%3 == stop_pos%3 && isStop(idx)  )
					{
					//there is another stop at this position in frame, not a good candidate
					return -1;
					}
				else if(!isATG(idx)) 
					{
					//nothing
					}
				else if(best==-1) {
					best= idx;
					}
				else if(best%3 != stop_pos%3) //not in frame
					{
					// nothing
					}
				else  {
					final KozakSequence kozak1= getKozakContext(best);
					final KozakSequence kozak2= getKozakContext(idx);
					if(kozak2.getStrength().compareTo(kozak1.getStrength())<0) {
						best= idx;
						}
					}
				}
			
			return best;
			}
		
		/** find best stop starting from ATG. 'best' is first in frame with atg */
		int findBestStop(final int ATG_pos) {
			if(isStop(ATG_pos)) throw new IllegalStateException("cannot be a stop !");
			int stop=-1;
			for(int idx=ATG_pos; idx+2 <  this.length();idx++) {
				if(isStop(idx)) {
					stop=idx;
					if(stop%3==ATG_pos%3) return stop;
					}
				}
			return stop;
			}
		
		protected boolean isPositiveStrand() {
			return this.getGene().isPositiveStrand();
			}
		@Override
		public String getContig() {
			return getGene().getContig();
			}
		@Override
		public int getStart() {
			return 1 + ( this.isPositiveStrand() ? this.getGene().getTxStart():this.getGene().getCdsEnd());
			}
		@Override
		public int getEnd() {
			return  ( this.isPositiveStrand() ? this.getGene().getCdsStart():this.getGene().getTxEnd());
			}
		/** return kozak context at ATG position . Return null if cannot build context **/
		KozakSequence getKozakContext(int atg0) {
			return new KozakSequence(this,atg0);
			}
		
		/** return peptide between pos_atg and post_stop or '.' if they're not in frame */
		String translate(int pos_ATG,int pos_stop) {
			if(pos_ATG>=pos_stop) throw new IllegalStateException("atg"+pos_ATG+">stop="+pos_stop);
			if(pos_ATG%3 != pos_stop%3) return ".";
			final StringBuilder sb=new StringBuilder();
			for(int i=pos_ATG;i<pos_stop;i+=3) {
				sb.append(VcfScanUpstreamOrf.this.geneticCode.translate(
						charAt(i+0),
						charAt(i+1),
						charAt(i+2))
						);
				}
			return sb.toString();
			}
		/** add annotation about stop */
		protected void annotStop(final StringBuilder sb,final int pos_ATG,final int pos_stop) {
			if(pos_stop==-1) {
				sb.append("stop:no");
				return;
				}
			
			if(pos_ATG!=-1)
				{
				sb. append("stop-frame:").
					append(pos_stop%3 == pos_ATG%3?"in-frame-stop":"not-in-frame-stop").
					append("|stop-pos:").
					append(1+this.convertOrfToGenomicIndex0(pos_stop)).
					append("|atg-stop:").
					append(pos_stop-pos_ATG).
					append("|pep:").
					append(pos_stop%3 == pos_ATG%3?translate(pos_ATG,pos_stop):".")
					;
				}
			}

		}
	
	private class UpstreamORF
		extends AbstractUTRSequence
		{
		final KnownGene knownGene; 
		private final int genomic_indexes0[];
		private final char base_buffer[];
		
		UpstreamORF(final KnownGene knownGene) {
			this.knownGene=knownGene;
			final List<Integer> L=new ArrayList<>();
			if(knownGene.isPositiveStrand())
				{				
				for(int exon_index=0;exon_index< knownGene.getExonCount();++exon_index) {
					if(knownGene.getExonStart(exon_index) >= knownGene.getCdsStart()) break;
					int beg = knownGene.getExonStart(exon_index);
					final int end = Math.min(knownGene.getExonEnd(exon_index),knownGene.getCdsStart());
					while(beg<end)
						{
						L.add(beg);
						++beg;
						}
					}
				}
			else
				{
				for(int exon_index=0;exon_index< knownGene.getExonCount();++exon_index) {
					if(knownGene.getExonEnd(exon_index) < knownGene.getCdsEnd()) continue;
					int beg = Math.max(knownGene.getExonStart(exon_index),knownGene.getCdsEnd());
					final int end = knownGene.getExonEnd(exon_index);
					while(beg<end)
						{
						L.add(beg);
						++beg;
						}
					}
				//Collections.reverse(L); no , will use binary search later 
				}
			this.genomic_indexes0 = L.stream().mapToInt(Integer::intValue).toArray();
			this.base_buffer = new char[this.genomic_indexes0.length];
			Arrays.fill(this.base_buffer, '\0');
			}
		
		
		@Override
		KnownGene getGene() {
			return this.knownGene;
			}
		
		@Override
		public final int length() {
			return this.genomic_indexes0.length;
			}
		
		boolean containsGenomicPos0(int genomic0) {
			final int x= Arrays.binarySearch(this.genomic_indexes0, genomic0);
			return x>=0 && x < this.genomic_indexes0.length;
			}
		
		int convertGenomicIndex0ToOrf(int g0) {
			final int x= Arrays.binarySearch(this.genomic_indexes0, g0);
			if( x<0 || x >= this.genomic_indexes0.length || this.genomic_indexes0[x]!=g0)
				{
				throw new IndexOutOfBoundsException("idx="+g0 +" length="+genomicSequence.length());
				}
			if(isPositiveStrand()) {
				return x;
				}
			else
				{
				return (this.genomic_indexes0.length-1)-x;
				}
			}
		
		@Override
		protected int convertOrfToGenomicIndex0(final int index) {
			if(index<0 || index>=this.length()) throw new IndexOutOfBoundsException("idx="+index +" length="+this.length());
			if(isPositiveStrand()) {
				return this.genomic_indexes0[index];
				}
			else
				{
				return this.genomic_indexes0[(this.genomic_indexes0.length-1)-index];
				}
			}

		@Override
		public final char charAt(final int index) {
			if(this.base_buffer[index]!='\0')
				{
				return this.base_buffer[index];
				}
			
			final int genomic_pos0 = this.convertOrfToGenomicIndex0(index);
			char c;
			if(isPositiveStrand()) {
				c= genomicSequence.charAt(genomic_pos0);
				}
			else {
				c= AcidNucleics.complement(genomicSequence.charAt(genomic_pos0));
				}
			c= Character.toUpperCase(c);
			this.base_buffer[index] = c;
			return c;
			}
				
		}
	
	private class MutatedUTR extends AbstractUTRSequence
		{
		final UpstreamORF delegate;
		final int genomic_index0_of_variant;
		final char alt_base;
		final int alt_position0;
		final Set<String> denovo_atg_set = new LinkedHashSet<>();
		final Set<String> remove_atg_set = new LinkedHashSet<>();
		final Set<String> remove_stop_set = new LinkedHashSet<>();
		final Set<String> denovo_stop_set = new LinkedHashSet<>();
		final Set<String> kozak_alterations_set = new LinkedHashSet<>();
		
		MutatedUTR(final UpstreamORF delegate,final VariantContext ctx,final int alt_index) {
			this.delegate=delegate;
			this.genomic_index0_of_variant = ctx.getStart()-1;
			char c=Character.toUpperCase((char)ctx.getAlleles().get(alt_index).getBases()[0]);
			this.alt_base= delegate.isPositiveStrand() ? c: AcidNucleics.complement(c);
			this.alt_position0 = delegate.convertGenomicIndex0ToOrf(this.genomic_index0_of_variant);
			if(delegate.convertOrfToGenomicIndex0(this.alt_position0)!=this.genomic_index0_of_variant) throw new IllegalStateException();
			}
		
		@Override
		protected int convertOrfToGenomicIndex0(final int index) {
			return this.delegate.convertOrfToGenomicIndex0(index);
			}
		
		@Override
		KnownGene getGene() {
			return this.delegate.getGene();
			}
		@Override
		public char charAt(int index) {
			return  index==this.alt_position0?
					this.alt_base:
					this.delegate.charAt(index);
			}
		
		@Override
		public int length() {
			return delegate.length();
			}
		
		boolean isChanging() {
			return !this.denovo_atg_set.isEmpty() ||
				!this.denovo_stop_set.isEmpty() ||
				!this.remove_atg_set.isEmpty() ||
				!this.remove_stop_set.isEmpty() ||
				!this.kozak_alterations_set.isEmpty()
				;
			}
		
		/* return atg shift position if alt base creates a new ATG: 0,1,2 or -1 if there is no new ATG */
		private StringBuilder getAnnotPrefix(final int pos_ATG) {
			final int dist_from_cap = pos_ATG;
			final int dist_to_cds_start = this.length() - pos_ATG ;

			
			return new StringBuilder("transcript:").
				append(this.getGene().getName()).
				append("|strand:").
				append(this.isPositiveStrand()?"+":"-").
				append("|utr-start:").
				append(this.getStart()).
				append("|utr-end:").
				append(this.getEnd()).
				append("|alt:").
				append(this.alt_base).
				append("|atg-pos:").
				append(1+this.delegate.convertOrfToGenomicIndex0(pos_ATG)).
				append("|cap-atg:").
				append(dist_from_cap).
				append("|atg-cds:").
				append(dist_to_cds_start).
				append("|atg-frame:").
				append( dist_to_cds_start % 3 !=0 ? "atg-out-of-cds-frame" :  "atg-in-cds-frame").
				append("|")
				;
			}

		
		void invoke() {
			if(!disable_atg_delete) removeATG();
			if(!disable_stop_delete) removeStop();
			if(!disable_atg_create) findDeNovoATG();
			if(!disable_stop_create) deNovoStop();
			if(!disable_kozak_alteration) kozakAlteration();
			}
		
		void removeStop() {
			for(int i=0;i<3;i++) {
				// stop removed
				if(!this.isStop(this.alt_position0 - i /* yes minus */)) continue;
				if(this.delegate.isStop(this.alt_position0 - i /* yes minus */)) continue;
				
				final int pos_stop = this.alt_position0 - i;
				
				final int pos_ATG = this.delegate.findBestATG(pos_stop);
				if(pos_ATG==-1) continue;
				if(pos_ATG>=pos_stop) throw new IllegalStateException(""+pos_ATG+"/"+pos_stop);
				
				
				
				final KozakSequence kozakSequence = this.getKozakContext(pos_ATG);
				if(acceptKozak(kozakSequence)) {
					continue;
					}
				final StringBuilder sb =getAnnotPrefix(pos_ATG);
				kozakSequence.annotate(sb);
				this.delegate.annotStop(sb, pos_ATG, pos_stop);
				
				this.remove_stop_set.add(sb.toString());
				}
			}
		
		void deNovoStop() {
			for(int i=0;i<3;i++) {
				// stop added
				if(!this.isStop(this.alt_position0 - i /* yes minus */)) continue;
			    if(this.delegate.isStop(this.alt_position0 - i /* yes minus */)) continue;
			    
				final int pos_stop = this.alt_position0 - i;
				
				final int pos_ATG = this.delegate.findBestATG(pos_stop);
				if(pos_ATG==-1) continue;
				final KozakSequence kozakSequence = this.getKozakContext(pos_ATG);
				if(!acceptKozak(kozakSequence)) continue;
				
				final StringBuilder sb =getAnnotPrefix(pos_ATG);
				
				kozakSequence.annotate(sb);
				
				this.annotStop(sb, pos_ATG, pos_stop);
				
				this.denovo_stop_set.add(sb.toString());
				}
			}

		
		
		void removeATG() {
			for(int i=0;i<3;i++) {
				if(!this.delegate.isATG(this.alt_position0 - i /* yes minus */)) continue;
				if(this.isATG(this.alt_position0 - i)) continue;
				
				
				final int pos_A = this.alt_position0 - i;
				final KozakSequence kozakSequence = this.delegate.getKozakContext(pos_A);
				if(!acceptKozak(kozakSequence)) continue;
				
				final int stop = this.findBestStop(pos_A);
				
				final StringBuilder sb =getAnnotPrefix(pos_A);
					
				kozakSequence.annotate(sb);
				annotStop(sb, pos_A, stop);
				
				this.remove_atg_set.add(sb.toString());
				}
			}
		
		
		void findDeNovoATG() {
			for(int i=0;i<3;i++) {
				//must be atg denovo
				if(!this.isATG(this.alt_position0 - i /* yes minus */)) continue;					
				if(this.delegate.isATG(this.alt_position0 - i /* yes minus */)) continue;					
				final int pos_A = this.alt_position0 - i;
				
				final KozakSequence kozakSequence = this.getKozakContext(pos_A);
				if(!acceptKozak(kozakSequence)) continue;
				final int stop = this.findBestStop(pos_A);
				
				final StringBuilder sb =getAnnotPrefix(pos_A);
				kozakSequence.annotate(sb);
				this.annotStop(sb,pos_A,stop);
				this.denovo_atg_set.add(sb.toString());
				}
			}
		
		void kozakAlteration() {
			for(int i=0;i< KOZAK_LENGTH;i++) {
				final int kozak_start = this.alt_position0 - i;
				if(kozak_start<0 || kozak_start+KOZAK_LENGTH>=this.length()) continue;
				final int pos_A = kozak_start+ KOZAK_ATG;
				final KozakSequence kozak1 = this.delegate.getKozakContext(pos_A);
				if(kozak1.getStrength()==KozakStrength.nil) continue;
				
				final KozakSequence kozak2 = this.getKozakContext(pos_A);
				if(kozak2.getStrength()==KozakStrength.nil) continue;
				
				if(!acceptKozak(kozak1) && !acceptKozak(kozak2)) {
					continue;
					}
				
				final StringBuilder sb = getAnnotPrefix(pos_A);
				
				sb.
					append("kozak-ref-seq:").
					append(kozak1.toString()).
					append("|kozak-ref-strength:").
					append(kozak1.getStrength().name()).
					append("|kozak-alt-seq:").
					append(kozak2.toString()).
					append("|kozak-alt-strength:").
					append(kozak2.getStrength().name()).
					append("|")
					;
				final int stop = this.findBestStop(pos_A);
				annotStop(sb, pos_A, stop);
				
				this.kozak_alterations_set.add(sb.toString());
				}
			}
		}
	
	private boolean acceptKozak(final KozakSequence k) {
		if(this.only_strong_kozak && !KozakStrength.Strong.equals(k.getStrength())) return false;
		return true;
	}
	
	/** rconvert transcript to interval return null if there is no UTR */
	private Interval kgToUTRInterval(final KnownGene kg) {
		if(kg.isPositiveStrand()) {
			//no utr
			if(kg.getStart()>=kg.getCdsStart()) return null;
			return new Interval(kg.getContig(),kg.getStart()+1,kg.getCdsStart(),false,kg.getName());
			}
		else  if(kg.isNegativeStrand())
			{
			//no utr
			if(kg.getCdsEnd()>=kg.getEnd()) return null;
			return new Interval(kg.getContig(),kg.getCdsEnd()/* no +1 */,kg.getEnd(),true,kg.getName());
			}
		else
			{
			LOG.info("bad strand in "+kg);
			return null;
			}
		}
		
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iter,
		final VariantContextWriter out
		) {
		try {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);
						
			
			LOG.info("Loading "+this.knownGeneUri);
			try(BufferedReader br= IOUtils.openURIForBufferedReading(this.knownGeneUri)) {
				String line;
				final Set<String> unknownContigs = new HashSet<>();
				final CharSplitter tab=CharSplitter.TAB;
				while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line))continue;
					
					final String tokens[]=tab.split(line);
					final KnownGene kg=new KnownGene(tokens);
					if(this.plus_strand_only && !kg.isPositiveStrand()) continue;
					
					final String refContig=this.refCtgNameConverter.apply(kg.getContig());
					if(StringUtils.isBlank(refContig)) {
						if(unknownContigs.add(kg.getContig())) {
							LOG.warn("unknown contig: "+kg.getContig());
							}
						continue;
						}
					kg.setChrom(refContig);
					
					
					final Interval interval = kgToUTRInterval(kg);
					if(interval==null) continue;
					if(kg.isPositiveStrand()) {}
					List<KnownGene> L =  this.knownGeneMap.get(interval);
					if(L==null) {
						L=new ArrayList<>();
						this.knownGeneMap.put(interval,L);
						}
					if(this.canonical_utr) {
						if(L.stream().
							map(K->kgToUTRInterval(K)).
							anyMatch(R->R.isNegativeStrand()==interval.isNegativeStrand() && R.equals(interval))
							) continue;
						}
					L.add(kg);
					}
				
				LOG.info("number of transcripts :"+this.knownGeneMap.values().stream().flatMap(L->L.stream()).count());
				}
  
			if(this.knownGeneMap.isEmpty()) {
				LOG.error("no transcripts found in "+this.knownGeneUri);
				return -1;
				}

			/** build vcf header */						
			final VCFHeader header= iter.getHeader();
			header.addMetaDataLine(this.infoAddATG);
			header.addMetaDataLine(this.infoAddSTOP);
			header.addMetaDataLine(this.infoDelATG);
			header.addMetaDataLine(this.infoDelSTOP);
			header.addMetaDataLine(this.infoKozakAlteration);
			JVarkitVersion.getInstance().addMetaData(this, header);
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			
			out.writeHeader(header);
			
			while(iter.hasNext()) {
				final VariantContext ctx = progress.apply(iter.next());
				if(!ctx.isVariant()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				if(!ctx.isSNP()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				final String refContig = this.refCtgNameConverter.apply(ctx.getContig());
				if(StringUtils.isBlank(refContig)) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				
				/* new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig)) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, refContig);
					}

				
				final Interval interval = new Interval(refContig,ctx.getStart(),ctx.getEnd()); 
				final List<KnownGene> kgGenes = this.knownGeneMap.getOverlapping(interval).
						stream().
						flatMap(C->C.stream()).
						collect(Collectors.toList());
				
				if(kgGenes.isEmpty()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				
				final List<UpstreamORF> uorfs = kgGenes.
						stream().
						map(KG->new UpstreamORF(KG)).
						filter(KG->KG.containsGenomicPos0(ctx.getStart()-1)).
						sorted((A,B)->Integer.compare(A.knownGene.getStart(),B.knownGene.getStart())).
						collect(Collectors.toList());
				
				if(uorfs.isEmpty()) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
								
				final List<MutatedUTR> mutatedList = new ArrayList<>();
				for(final UpstreamORF uorf: uorfs) {
					for(int alt_idx=1 /* 0==REF */;alt_idx< ctx.getAlleles().size();++alt_idx) {
						final Allele alt_allele = ctx.getAlleles().get(alt_idx);
						if(alt_allele.isSymbolic() || !alt_allele.isCalled() || alt_allele.length()!=1) continue;
						if(!AcidNucleics.isATGC(alt_allele.getDisplayBases()[0])) continue;
						final MutatedUTR mutated = new MutatedUTR(uorf, ctx, alt_idx);
						mutated.invoke();
						mutatedList.add(mutated);
						}
					}
				
				if(mutatedList.isEmpty() || mutatedList.stream().noneMatch(M->M.isChanging())) {
					if(!this.print_uorf_only) out.add(ctx);
					continue;
					}
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				
				List<String> ann = new ArrayList<>(
							mutatedList.
							stream().
							flatMap(M->M.remove_atg_set.stream()).
							collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoDelATG.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.remove_stop_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoDelSTOP.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.denovo_atg_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoAddATG.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.denovo_stop_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoAddSTOP.getID(), ann);
					}
				
				ann = new ArrayList<>(
						mutatedList.
						stream().
						flatMap(M->M.kozak_alterations_set.stream()).
						collect(Collectors.toCollection(LinkedHashSet::new)));
				if(!ann.isEmpty()) {
					vcb.attribute(this.infoKozakAlteration.getID(), ann);
					}
				
				out.add(vcb.make());
				}
			
			progress.close();
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		return super.doVcfToVcf(args, this.outputFile);
		}
	
	public static void main(final String[] args) {
		new VcfScanUpstreamOrf().instanceMainWithExit(args);
	}
}