/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;

import com.github.lindenb.jvarkit.bed.BedCoordMath;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;

import htsjdk.tribble.annotation.Strand;

/**
 * Implementation of UcscTranscript in the genepred format
 * @author lindenb
 *
 */
class UcscTranscriptImpl implements UcscTranscript {
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
	
	/** an abstract RNA : start() and end() on genomic sequence need to be defined*/
	abstract class RNAImpl extends DelegateCharSequence implements UcscTranscript.RNA
		{
		/** cache length */
		private Integer _length=null;
		RNAImpl(final CharSequence sequence)
			{
			super(sequence);
			}

		/** get Gene associated to this RNA */
		@Override
		public final UcscTranscriptImpl getTranscript()
			{
			return UcscTranscriptImpl.this;
			}
		/** start of mRNA (could be transcription or traduction */
		protected abstract int start();
		/** end of mRNA (could be transcription or traduction */
		protected abstract int end();

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

		/** convert the genomic position to the position in the RNA, return -1 if RNA not in genomic pos */
		OptionalInt convertGenomic0ToRnaCoordinate0(int genomicPos0)
			{
			int rnaPos0=0;
			if(getTranscript().isPositiveStrand())
				{
				for(final Exon ex:getTranscript().getExons())
					{
					if(this.start()>=ex.getEnd()) continue;
					if(this.end()<=ex.getStart()) break;
					final int beg=Math.max(this.start(), ex.getStart());
					final int end=Math.min(this.end(), ex.getEnd());
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
					if(this.start()>=ex.getEnd()) break;
					if(this.end()<=ex.getStart()) continue;
					final int beg=Math.max(this.start(), ex.getStart());
					final int end=Math.min(this.end(), ex.getEnd());
					if(beg<= genomicPos0 && genomicPos0<end)
						{
						return OptionalInt.of(((end-1)-genomicPos0) + rnaPos0);
						}
					rnaPos0+=(end-beg);
					}
				}
			return OptionalInt.empty();
			}
		
		int convertToGenomicCoordinate(int rnaPos0)
			{
			if(rnaPos0<0) throw new IllegalArgumentException("negative index:"+rnaPos0);
			if(rnaPos0>=this.length()) throw new IndexOutOfBoundsException("out of bound index:"+rnaPos0+"<"+this.length());
			
			if(isPositiveStrand())
				{
				for(Exon ex:getTranscript().getExons())
					{
					if(this.start()>=ex.getBedEnd()) continue;
					if(this.end()<=ex.getBedStart()) break;
					final int beg=Math.max(this.start(), ex.getBedStart());
					final int len= BedCoordMath.getOverlap(this.start(), this.end(), ex.getBedStart(), ex.getBedEnd());
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
					if(this.start() >= ex.getBedEnd()) break;
					if(this.end() <= ex.getBedStart()) continue;
					final int end=Math.min(this.end(), ex.getBedEnd());
					final int len= BedCoordMath.getOverlap(this.start(), this.end(), ex.getBedStart(), ex.getBedEnd());
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
					if(this.start() >= ex.getBedEnd()) continue;
					if(this.end() <= ex.getBedStart()) break;
					this._length+= BedCoordMath.getOverlap(
							this.start(),this.end(),
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
			int n= convertToGenomicCoordinate(index0);
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
		protected  final int start()
			{
			return getTranscript().getTxStart() -1;
			}
		@Override
		protected  final int end()
			{
			return getTranscript().getTxEnd();
			}
		
		public CodingRNAImpl getCodingRNA() {
			if(buffer != null) return buffer;
			if(!hasCodingRNA()) throw new IllegalStateException("mRNA does not code for a protein");
			buffer = new CodingRNAImpl(this);
			return buffer;
			}
		
		public List<CodingRNAImpl> getUpstreamORFs() {
			if(!hasCodingRNA()) throw new IllegalStateException("mRNA does not code for a protein");
			int rnaX,rnaY;
			if(isPositiveStrand()) {
				rnaX = 0;
				rnaY = convertGenomic0ToRnaCoordinate0(getCdsStart()-1).getAsInt();
				}
			else
				{
				rnaX = convertGenomic0ToRnaCoordinate0(getCdsEnd()).getAsInt();
				rnaY = length();
				}
			return getORFs(rnaX,rnaY);
			}
		
		protected List<CodingRNAImpl> getORFs(final int rnaX1, final int rnaX2) {
			List <CodingRNAImpl> L = null;
			int i = rnaX1;
			while(i+2 < length() && i+2 < rnaX2) {
				if(isATG(i)) {
					int atg0 = i;
					i+=3;
					while(i+2 < rnaX2 && i+2< length() && !isStop(i)) {
						
						i+=3;
						}
					if(L==null) L = new ArrayList<>();
					int a = convertToGenomicCoordinate(atg0);
					int b = convertToGenomicCoordinate(i);
					if(isNegativeStrand()) {
						int c = a;
						a = b-2;
						b = c-2;
						}
					L.add(new CodingRNAImpl(this,a,b));
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
		protected  final int start()
			{
			return this._start0;
			}
		@Override
		protected  final int end()
			{
			return this._end0;
			}
		int convertToRNAMessengerCoordinate(int cdna0) {
			return -1;//TODO
			}
		@Override
		public KozakSequence getKozakSequence() {
			final int atg0 = convertToRNAMessengerCoordinate(0);
			return new KozakSequence(getMessengerRNA().getDelegate(),atg0);
			}
		
		public Peptide getPeptide()
			{
			return new PeptideImpl(GeneticCode.getStandard(), this);
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
				final int idx1=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+0);
				final int idx2=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+1);
				final int idx3=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+2);
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
		public int[] convertToGenomicCoordinates(int pepPos0)
			{
			if(pepPos0<0) throw new IndexOutOfBoundsException("negative offset : "+pepPos0);
			if(pepPos0>=this.length()) throw new IndexOutOfBoundsException(" idx="+pepPos0+" and length="+this.length());
			
			return new int[]{
				getCodingRNA().convertToGenomicCoordinate(pepPos0*3+0),	
				getCodingRNA().convertToGenomicCoordinate(pepPos0*3+1),	
				getCodingRNA().convertToGenomicCoordinate(pepPos0*3+2)
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
	
	
	UcscTranscriptImpl() {
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
	public MessengerRNA getMessengerRNA(CharSequence chromosomeSequence) {
		return new MessengerRNAImpl(chromosomeSequence);
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
			sb.append(" cdsStart:").append(getCdsStart()).
			append(" cdsEnd:").append(getCdsEnd());
			}
			
		return sb.toString();
		}
	}
