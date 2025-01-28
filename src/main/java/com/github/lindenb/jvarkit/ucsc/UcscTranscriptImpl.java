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
package com.github.lindenb.jvarkit.ucsc;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.bed.BedCoordMath;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;

/**
 * Implementation of UcscTranscript in the genepred format
 * @author lindenb
 *
 */
class UcscTranscriptImpl implements UcscTranscript {
	private static final boolean DEBUG = Boolean.parseBoolean(System.getProperty("ucsc.transcript.debug","false"));
	private static final Logger LOG = Logger.build(UcscTranscriptImpl.class).make();
	
	String contig;
	int txStart;
	int txEnd;
	int cdsStart;
	int cdsEnd;
	char strand;
	int[] exonStarts;
	int[] exonEnds;
	String name;
	String name2 = null;
	OptionalInt score = OptionalInt.empty();
	Map<String,String> metadata = Collections.emptyMap();
	
		
	
	
	/** an abstract RNA : start() and end() on genomic sequence need to be defined*/
	abstract class RNAImpl extends DelegateCharSequence implements UcscTranscript.RNA
		{
		/** cache length */
		private Integer _length=null;
		RNAImpl(final CharSequence sequence)
			{
			super(sequence);
			}
		
		@Override
		public Strand getStrand() {
			return getTranscript().getStrand();
			}


		@Override
		public int getStart() {
			return  getBedStart()+1;
			}
		@Override
		public int getEnd() {
			return getBedEnd();
			}
		
		/** get Gene associated to this RNA */
		@Override
		public /** final no : overriden by internal ORF */ UcscTranscriptImpl getTranscript()
			{
			return UcscTranscriptImpl.this;
			}

		private boolean hasCodon(final int pos0, char b1,char b2, char b3) {
			if(pos0 < 0 || pos0 +2 >= this.length()) return false;
			char c = Character.toUpperCase(charAt(pos0  ));
			if(c!=b1) return false;
			c = Character.toUpperCase(charAt(pos0+1));
			if(c!=b2) return false;
			c = Character.toUpperCase(charAt(pos0+2));
			if(c!=b3) return false;
			return true;
		}

		protected boolean isATG(final int pos0) {
			return hasCodon(pos0,'A','T','G');
			}
		protected boolean isStop(final int pos0) {
			return hasCodon(pos0,'T','A','G') ||
				hasCodon(pos0,'T','A','A') ||
				hasCodon(pos0,'T','G','A')
				;
			}

		void checkCoordinates() {
			if(!DEBUG) return;
			for(int i=0;i< length();i++) {
				int g0 = convertToGenomic0Coordinate0(i);
				OptionalInt r0 = convertGenomic0ToRnaCoordinate0(g0);
				if(!r0.isPresent()|| r0.getAsInt()!=i) throw new IllegalStateException();
				}
			}


		/** convert the genomic position to the position in the RNA, return -1 if RNA not in genomic pos */
		@Override
		public OptionalInt convertGenomic0ToRnaCoordinate0(int genomicPos0)
			{
			if(DEBUG) LOG.debug("convertGenomic0ToRnaCoordinate0 g0="+genomicPos0+" for "+ getTranscript());
			int rnaPos0=0;
			if(getTranscript().isPositiveStrand())
				{
				for(final Exon ex:getTranscript().getExons())
					{
					if(this.getBedStart() >= ex.getBedEnd()) continue;
					if(this.getBedEnd() <= ex.getBedStart()) break;
					final int beg = Math.max(this.getBedStart(), ex.getBedStart());
					final int end = Math.min(this.getBedEnd(), ex.getBedEnd());
					if(beg<= genomicPos0 && genomicPos0<end)
						{
						return OptionalInt.of((genomicPos0-beg)+rnaPos0);
						}
					rnaPos0+=(end-beg);
					}
				}
			else
				{
				for(int idx=getTranscript().getExonCount()-1;idx>=0;idx--)
					{
					final Exon ex=getExon(idx);
					if(DEBUG) System.err.println(ex.getName()+" ex.start:"+ex.getBedStart()+"-"+ex.getBedEnd()+" "+rnaPos0+" "+this.getBedStart()+"-"+this.getBedEnd());
					if(this.getBedStart()>=ex.getBedEnd()) break;
					if(this.getBedEnd()<=ex.getBedStart()) continue;
					final int beg = Math.max(this.getBedStart(), ex.getBedStart());
					final int end = Math.min(this.getBedEnd(), ex.getBedEnd());
					if(beg<= genomicPos0 && genomicPos0<end)
						{
						return OptionalInt.of(((end-1)-genomicPos0) + rnaPos0);
						}
					rnaPos0+=(end-beg);
					}
				}
			return OptionalInt.empty();
			}
		@Override
		public int convertToGenomic0Coordinate0(int rnaPos0)
			{
			if(rnaPos0<0) throw new IllegalArgumentException("negative index:"+rnaPos0);
			if(rnaPos0>=this.length()) throw new IndexOutOfBoundsException("out of bound index:"+rnaPos0+"<"+this.length());
			
			if(isPositiveStrand())
				{
				for(Exon ex:getTranscript().getExons())
					{
					if(this.getBedStart()>=ex.getBedEnd()) continue;
					if(this.getBedEnd()<=ex.getBedStart()) break;
					final int beg=Math.max(this.getBedStart(), ex.getBedStart());
					final int len= BedCoordMath.getOverlap(this.getBedStart(), this.getBedEnd(), ex.getBedStart(), ex.getBedEnd());
					if(rnaPos0<len)
						{
						return beg+rnaPos0;
						}
					rnaPos0-=len;
					}
				}
			else
				{
				for(int idx=getTranscript().getExonCount()-1;idx>=0;idx--)
					{
					final Exon ex=getExon(idx);
					if(this.getBedStart() >= ex.getBedEnd()) break;
					if(this.getBedEnd() <= ex.getBedStart()) continue;
					final int end=Math.min(this.getBedEnd(), ex.getBedEnd());
					final int len= BedCoordMath.getOverlap(this.getBedStart(), this.getBedEnd(), ex.getBedStart(), ex.getBedEnd());
					if(rnaPos0<len)
						{
						return (end-1)-rnaPos0;
						}
					rnaPos0-=len;
					}
				}
			throw new IllegalArgumentException("Illegal State");
			}
		
		
		
		@Override
		public int length()
			{
			if(this._length==null)
				{
				this._length=0;
				for(final Exon ex:getTranscript().getExons())
					{
					if(this.getBedStart() >= ex.getBedEnd()) continue;
					if(this.getBedEnd() <= ex.getBedStart()) break;
					this._length+= BedCoordMath.getOverlap(
							this.getBedStart(),this.getBedEnd(),
							ex.getBedStart(),ex.getBedEnd()
							);
					}
				}
			return this._length;
			}
		@Override
		public char charAt(int index0)
			{
			if(index0<0) throw new IllegalArgumentException("negative index:"+index0);
			if(index0>=this.length()) throw new IndexOutOfBoundsException("index:"+index0 +" < "+this.length());
			int n= convertToGenomic0Coordinate0(index0);
			if(n==-1) 	throw new IndexOutOfBoundsException("0<=index:="+index0+"<"+length());
			if(isPositiveStrand())
				{
				return getDelegate().charAt(n);
				}
			else
				{
				return AcidNucleics.complement(getDelegate().charAt(n));
				}
			}
		}

	public class MessengerRNAImpl extends RNAImpl implements UcscTranscript.MessengerRNA {
		private CodingRNAImpl buffer = null;
		MessengerRNAImpl(final CharSequence sequence)
			{
			super(sequence);
			}
	
		@Override
		public final int getBedStart()
			{
			return getTranscript().getTxStart() -1;
			}
		@Override
		public final int getBedEnd()
			{
			return getTranscript().getTxEnd();
			}
		
		@Override
		public CodingRNAImpl getCodingRNA() {
			if(buffer != null) return buffer;
			if(!hasCodingRNA()) throw new IllegalStateException("mRNA does not code for a protein");
			buffer = new CodingRNAImpl(this);
			buffer.checkCoordinates();
			return buffer;
			}
		

		
		@Override
		public UntranslatedRNAImpl getUpstreamUntranslatedRNA() {
			if(!hasUpstreamUntranslatedRNA()) throw new IllegalStateException("No upstream UTR in "+ getTranscript().getTranscriptId());
			if(isPositiveStrand()) {
				return new UntranslatedRNAImpl(this, UcscTranscriptImpl.this.txStart,UcscTranscriptImpl.this.cdsStart);
				}
			else
				{
				return new UntranslatedRNAImpl(this, UcscTranscriptImpl.this.cdsEnd,UcscTranscriptImpl.this.txEnd);
				}
			}
		
		
		
		
		@Override
		public UntranslatedRNAImpl getDownstreamUntranslatedRNA() {
			if(!hasDownstreamUntranslatedRNA()) throw new IllegalStateException("No downstream UTR in "+ getTranscript().getTranscriptId());
			if(isPositiveStrand()) {
				return new UntranslatedRNAImpl(this, UcscTranscriptImpl.this.cdsEnd,UcscTranscriptImpl.this.txEnd);
				}
			else
				{
				return new UntranslatedRNAImpl(this, UcscTranscriptImpl.this.txStart,UcscTranscriptImpl.this.cdsStart);
				}
			}
		
		
		protected List<CodingRNA> getORFs(final int rnaX1, final int rnaX2) {
			List<CodingRNA> L = null;
			for(int phase=0;phase<3;++phase) {
				int i = rnaX1 + phase;
				while(i+2 < length() && i+2 < rnaX2) {
					if(isATG(i)) {
						int atg0 = i;
						i+=3;
						while(i+2< length() && i+2 < rnaX2) {
							if(isStop(i)) break;
							i+=3;
							}
						if(L==null) L = new ArrayList<>();
						final int a;
						final int b;
						if(isPositiveStrand()) {
							a = convertToGenomic0Coordinate0(atg0);
							b = convertToGenomic0Coordinate0(i);
							}
						else
							{
							a = convertToGenomic0Coordinate0(i);
							b = convertToGenomic0Coordinate0(atg0)+1;//one base past ATG
							}
						if(a>=b) throw new IllegalStateException("boum");
						final CodingRNAImpl coding = new MicroCodingRNAImpl(this,a,b);
						if(!coding.isATG(0)) throw new IllegalStateException("not an atg in "+coding+" "+getContig()+":"+a+"-"+b+" "+getStrand());
						L.add(coding);
						}
					else
						{
						i+=3;
						}
					}
				}
			return L==null?Collections.emptyList():L;
			}
		}


	/** same as RNA but without the UTR, starts is cdsStart, end is cdsEnd */
	public class CodingRNAImpl  extends RNAImpl implements UcscTranscript.CodingRNA
		{
		private final MessengerRNAImpl mRNA;
		private final int _start0;
		private final int _end0;
		CodingRNAImpl(final MessengerRNAImpl mRNA, int _start0,int _end0) {
			super(mRNA.getDelegate());
			this.mRNA = mRNA;
			this._start0 = _start0;
			this._end0 = _end0;
			}
		CodingRNAImpl(final MessengerRNAImpl mRNA) {
			this(mRNA, mRNA.getTranscript().getCdsStart()-1,mRNA.getTranscript().getCdsEnd());
			}
		@Override
		public MessengerRNAImpl getMessengerRNA() {
			return this.mRNA;
			}
		@Override
		public  final int getBedStart()
			{
			return this._start0;
			}
		@Override
		public  final int getBedEnd()
			{
			return this._end0;
			}
		@Override
		public int convertCoding0ToMessenger0(int cdna0) {
			final int g0 = convertToGenomic0Coordinate0(cdna0);
			return getMessengerRNA().convertGenomic0ToRnaCoordinate0(g0).getAsInt();
			}
		@Override
		public KozakSequence getKozakSequence() {
			final int atg0 = convertCoding0ToMessenger0(0);
			return new KozakSequence(getMessengerRNA(),atg0);
			}
		
		public Peptide getPeptide()
			{
			return new PeptideImpl(GeneticCode.getStandard(), this);
			}
		}


	class UntranslatedRNAImpl extends RNAImpl implements UcscTranscript.UntranslatedRNA {
		private final MessengerRNAImpl mRNA;
		private final int utrStart;
		private final int utrEnd;
		UntranslatedRNAImpl(final MessengerRNAImpl mRNA,int utrStart,int utrEnd) {
			super(mRNA.getDelegate());
			this.mRNA = mRNA;
			this.utrStart = utrStart;
			this.utrEnd = utrEnd;
			if(this.utrStart>=this.utrEnd) throw new IllegalArgumentException("utrStart>=utrEnd");
			}
		@Override
		public MessengerRNAImpl getMessengerRNA() {
			return this.mRNA;
			}
		@Override
		public  final int getBedStart()
			{
			return this.utrStart;
			}
		@Override
		public  final int getBedEnd()
			{
			return this.utrEnd;
			}
		
		
		
		public List<CodingRNA> getORFs() {			
			final int atg_genomic0 = getTranscript().getGenome0ATG();
			final OptionalInt atg0 = getMessengerRNA().convertGenomic0ToRnaCoordinate0(atg_genomic0);
			if(!atg0.isPresent()) {
				throw new IllegalStateException("Cannot get atg0 for "+getTranscript()+" g0="+atg_genomic0+" "+getMessengerRNA().hasUpstreamUntranslatedRNA());
				}
			return getMessengerRNA().getORFs(0, atg0.getAsInt());
			}
		}

	public class MicroCodingRNAImpl extends CodingRNAImpl {
		final UcscTranscriptImpl transcript2;
		MicroCodingRNAImpl(final MessengerRNAImpl mRNA, int _start0,int _end0) {
			super(mRNA,_start0,_end0);
			this.transcript2 = new UcscTranscriptImpl(mRNA.getTranscript()) {
				@Override
				public MessengerRNA getMessengerRNA(CharSequence chromosomeSequence) {
					return mRNA;
					}
				};
			this.transcript2.name += ".uORF."+(_start0+1)+"-"+(_end0);
			this.transcript2.cdsStart = _start0;
			this.transcript2.cdsEnd = _end0;
			}
		@Override
		public final UcscTranscriptImpl getTranscript() {
			return transcript2;
			}
		}
	
	/** implementation of UcscTranscript.Peptide */
	public class PeptideImpl extends DelegateCharSequence  implements UcscTranscript.Peptide
		{
		private Integer _length=null;
		private final GeneticCode geneticCode;
		PeptideImpl(final GeneticCode gc, final CodingRNA rna)
			{
			super(rna);
			this.geneticCode=gc;
			}
		@Override
		public CodingRNAImpl getCodingRNA()
			{
			return CodingRNAImpl.class.cast(getDelegate());
			}
		@Override
		public int length()
			{
			if(_length==null)
				{
				_length = getDelegate().length()/3;
				if(_length==0)
					{
					return 0;
					}
				final int idx1=getCodingRNA().convertToGenomic0Coordinate0((_length-1)*3+0);
				final int idx2=getCodingRNA().convertToGenomic0Coordinate0((_length-1)*3+1);
				final int idx3=getCodingRNA().convertToGenomic0Coordinate0((_length-1)*3+2);
				if(idx1==-1 || idx2==-1 || idx3==-1)
					{
					System.err.println("Bizarre pour "+getCodingRNA().getTranscript().getTranscriptId()+" "+getStrand());
					return _length;
					}
				
				final char last=this.getGeneticCode().translate(
						getCodingRNA().charAt((_length-1)*3+0),
						getCodingRNA().charAt((_length-1)*3+1),
						getCodingRNA().charAt((_length-1)*3+2)
						);
				if(!Character.isLetter(last)) _length--;
				}
			return _length;
			}
		
		public char[] getCodon(int pepPos0)
			{
			return new char[]{
				getCodingRNA().charAt(pepPos0*3+0),
				getCodingRNA().charAt(pepPos0*3+1),
				getCodingRNA().charAt(pepPos0*3+2)
				};
			}

		/** convert the genomic position to the position in the peptide, return -1 if peptide not in genomic pos */
		public OptionalInt convertGenomicToPeptideCoordinate(int genomicPos0)
			{
			final OptionalInt rnaIdx=getCodingRNA().convertGenomic0ToRnaCoordinate0(genomicPos0);
			if(!rnaIdx.isPresent()) return OptionalInt.empty();
			return OptionalInt.of(rnaIdx.getAsInt()/3);
			}

		@Override
		public int[] convertToGenomic0Coordinates(int pepPos0)
			{
			if(pepPos0<0) throw new IndexOutOfBoundsException("negative offset : "+pepPos0);
			if(pepPos0>=this.length()) throw new IndexOutOfBoundsException(" idx="+pepPos0+" and length="+this.length());
			
			return new int[]{
				getCodingRNA().convertToGenomic0Coordinate0(pepPos0*3+0),	
				getCodingRNA().convertToGenomic0Coordinate0(pepPos0*3+1),	
				getCodingRNA().convertToGenomic0Coordinate0(pepPos0*3+2)
				};
			}
		
		@Override
		public char charAt(final int pepPos0) {
			return this.getGeneticCode().translate(
					getCodingRNA().charAt(pepPos0*3+0),	
					getCodingRNA().charAt(pepPos0*3+1),	
					getCodingRNA().charAt(pepPos0*3+2)
					);
			}
		
		protected KozakSequence getKozakSequence(int mrna0) {
			return new KozakSequence(this,mrna0);
			}
		
		public GeneticCode getGeneticCode() {
			return this.geneticCode;
			}
		}


	private class UTRImpl extends UcscTranscript.UTR {
		private int start0;
		private int end0;
		UTRImpl(final UcscTranscript.Exon exon,int start0,int end0) {
			super(exon);
			this.start0 = start0;
			this.end0 = end0;
			if(this.start0>=this.end0) throw new IllegalArgumentException("start0<end0");
			}
		@Override
		public int getStart() {
			return start0+1;
			}
		@Override
		public int getEnd() {
			return end0;
			}
		@Override
		public boolean isUTR5() {
			if(isPositiveStrand()) {
				return BedCoordMath.overlaps(start0,end0,getTranscript().getTxStart()-1,getTranscript().getCdsStart()-1);
				}
			else
				{
				return BedCoordMath.overlaps(start0,end0,getTranscript().getCdsEnd(),getTranscript().getTxEnd());
				}
			}
		@Override
		public boolean isUTR3() {
			return !isUTR5();
			}
		@Override
		public String getName() {
			return isUTR5()?"UTR5":"UTR3";
			}
		}
	
	@Override
	public List<UTR> getUTRs() {
		if(!isProteinCoding()) return Collections.emptyList();
		if(this.txStart==this.cdsStart && this.txEnd==this.cdsEnd) return Collections.emptyList();
		final List<UTR> L = new ArrayList<>();
		for(int i=0;i< getExonCount();i++) {
			final Exon ex = getExon(i);
			if(BedCoordMath.overlaps(ex.getBedStart(),ex.getBedEnd(),this.txStart,this.cdsStart)) {
				L.add(new UTRImpl(ex,
					Math.max(ex.getBedStart(),this.txStart),
					Math.min(ex.getBedEnd(),this.cdsStart)
					));
				}
			if(BedCoordMath.overlaps(ex.getBedStart(),ex.getBedEnd(),this.cdsEnd,this.txEnd)) {
				L.add(new UTRImpl(ex,
					Math.max(ex.getBedStart(),this.cdsEnd),
					Math.min(ex.getBedEnd(),this.txEnd)
					));
				}
			}
		return L;
		}
	
	
	private class CodonImpl extends UcscTranscript.Codon {
		private int start1;
		private int end1;
		private boolean is_start_codon;
		CodonImpl(final Exon exon,int start1,int end1,boolean is_start_codon) {
			super(exon);
			this.start1 = start1;
			this.end1 = end1;
			this.is_start_codon=is_start_codon;
			}
		
		@Override
		public int getStart() {
			return start1;
			}
		@Override
		public int getEnd() {
			return end1;
			}

		@Override
		public boolean isStartCodon() {
			return is_start_codon;
		}

		@Override
		public String getName() {
			return isStartCodon()?"start_codon":"end_codon";
			}
		}
	
	/** get last 3 positions for a codon start/stop */
	private static  List<Integer> last3(IntStream stream) {
	    Deque<Integer> result = new ArrayDeque<>(3);
	    stream.forEachOrdered(x -> {
	        if (result.size() == 3) {
	            result.pop();
	        	}
	        result.add(x);
	    	});
	    return new ArrayList<>(result);
		}
	
	/** convert 3 successive position to codon */
	private static  List<Locatable> position_to_codons(String contig,List<Integer> positions) {
		if(positions.size()!=3) return Collections.emptyList();
		List<Locatable> L=new ArrayList<>();
		int p1 = positions.get(0);
		int p2 = positions.get(1);
		int p3 = positions.get(2);
		if(p1+1==p2 && p2+1==p3) {
			L.add(new SimpleInterval(contig,p1,p3));
			}
		else if(p1+1!=p2 && p2+1==p3) {
			L.add(new SimpleInterval(contig,p1,p1));
			L.add(new SimpleInterval(contig,p2,p3));
			}
		else if(p1+1==p2 && p2+1!=p3) {
			L.add(new SimpleInterval(contig,p1,p2));
			L.add(new SimpleInterval(contig,p3,p3));
			}
		else
			{
			L.add(new SimpleInterval(contig,p1,p1));
			L.add(new SimpleInterval(contig,p2,p2));
			L.add(new SimpleInterval(contig,p3,p3));
			}
		return L;
		}
	
	@Override
	public List<Codon> getCodons() {
		if(!isProteinCoding()) return Collections.emptyList();
		if(this.txStart==this.cdsStart && this.txEnd==this.cdsEnd) return Collections.emptyList();
		final List<CDS> cdss = getCDS();
		if(cdss.isEmpty()) return Collections.emptyList();
		
		final List<Codon> return_codon = new ArrayList<>();
		
		final List<Locatable> list5 = position_to_codons(getContig(), cdss.stream().
				flatMapToInt(CDS->IntStream.rangeClosed(CDS.getStart(), CDS.getEnd()).limit(3L)).
				limit(3L).
				boxed().
				collect(Collectors.toList())
				);
			
		final List<Locatable> list3 = position_to_codons(getContig(),last3(cdss.stream().
					flatMapToInt(CDS->IntStream.rangeClosed(CDS.getStart(), CDS.getEnd()))
					));
		
		for(int side=0;side<2;++side) {
			for(Locatable codonloc: (side==0?list5:list3)) {
				final Exon exon = getExons().stream().filter(EX->EX.overlaps(codonloc)).findFirst().get();
				final boolean is_start_codon = ((side==0&&this.isPositiveStrand()) || (side==1 && this.isNegativeStrand()));
				final CodonImpl codon = new CodonImpl(exon, codonloc.getStart(), codonloc.getEnd(), is_start_codon);
				return_codon.add(codon);
			}
		}
		
		return return_codon;
		}
	
	
	UcscTranscriptImpl() {
		}



	UcscTranscriptImpl(final UcscTranscriptImpl cp) {
		this.contig = cp.contig;
		this.txStart = cp.txStart;
		this.txEnd = cp.txEnd;
		this.cdsStart = cp.cdsStart;
		this.cdsEnd = cp.cdsEnd;
		this.strand = cp.strand;
		this.exonStarts = cp.exonStarts;
		this.exonEnds = cp.exonEnds;
		this.name = cp.name;
		this.name2 = cp.name2;
		this.score = cp.score;
		this.metadata = cp.metadata;
	}
	
	
	@Override
	public Strand getStrand() {
		return Strand.decode(this.strand);
		}

	@Override
	public boolean isNegativeStrand() {
		return this.strand == '-';
		}
	@Override
	public boolean isPositiveStrand() {
		return this.strand == '+';
		}
	
	
	@Override
	public String getContig() {
		return this.contig;
		}
	
	
	@Override
	public final boolean isProteinCoding() {
		return this.cdsStart < this.cdsEnd;
		}
	
	@Override
	public int getCdsStart() {
		if(!isProteinCoding()) throw new IllegalStateException("Cannot ask cdsStart for non-protein coding. "+getTranscriptId());
		return this.cdsStart+1;
		}
	
	@Override
	public int getCdsEnd() {
		if(!isProteinCoding()) throw new IllegalStateException("Cannot ask cdsend for non-protein coding. "+getTranscriptId());
		return this.cdsEnd;
		}

	
	@Override
	public String getTranscriptId() {
		return name;
		}
	
	@Override
	public String getName2() {
		return name2;
		}
	
	@Override
	public int hashCode() {
		int i = getTranscriptId().hashCode();
		i= i*31 + getContig().hashCode();
		i= i*31 + Integer.hashCode(txStart);
		i= i*31 + Integer.hashCode(txEnd);
		return i;
		}
	
	@Override
	public int getBedStart() {
		return txStart;
		}
	@Override
	public int getBedEnd() {
		return txEnd;
		}
	
	@Override
	public int getTxStart() {
		return txStart+1;
		}
	
	@Override
	public int getTxEnd() {
		return txEnd;
		}
	
	@Override
	public int getExonCount() {
		return this.exonStarts.length;
		}
	
	@Override
	public int getExonStart(int idx) {
		return this.exonStarts[idx]+1;
		}
	
	@Override
	public int getExonEnd(int idx) {
		return this.exonEnds[idx];
		}
	
	@Override
	public MessengerRNA getMessengerRNA(final CharSequence chromosomeSequence) {
		final MessengerRNAImpl rna= new MessengerRNAImpl(chromosomeSequence);
		if(DEBUG) rna.checkCoordinates();
		return rna;
		}
	
	@Override
	public String toString() {
		final StringBuilder sb= new StringBuilder().
			append("contig:").append(getContig()).
			append(" name:").append(getTranscriptId()).
			append(" strand:").append(getStrand()).
			append(" txStart:").append(getTxStart()).
			append(" txEnd:").append(getTxEnd());
		if(isProteinCoding()) {
			sb.append(" coding:true");
			sb.append(" cdsStart:").append(getCdsStart()).
			append(" cdsEnd:").append(getCdsEnd());
			}
			
		return sb.toString();
		}
	}
