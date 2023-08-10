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

import java.util.OptionalInt;

import com.github.lindenb.jvarkit.bed.BedCoordMath;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;

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
	
	
	/** same as RNA but without the UTR, starts is cdsStart, end is cdsEnd */
	public class CodingRNAImpl  extends RNAImpl implements UcscTranscript.CodingRNA
		{
		CodingRNAImpl(final CharSequence sequence)
			{
			super(sequence);
			}
		@Override
		protected  final int start()
			{
			return getTranscript().getCdsStart() -1;	
			}
		@Override
		protected  final int end()
			{
			return getTranscript().getCdsEnd();	
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
		PeptideImpl(final GeneticCode gc,CodingRNA rna)
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
		
		public GeneticCode getGeneticCode() {
			return this.geneticCode;
			}
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
	public boolean isProteinCoding() {
		return this.cdsStart < this.cdsEnd;
		}
	
	@Override
	public int getCdsStart() {
		if(!isProteinCoding()) throw new IllegalStateException("Cannot ask cdsStart for non-protein coding");
		return this.cdsStart+1;
		}
	
	@Override
	public int getCdsEnd() {
		if(!isProteinCoding()) throw new IllegalStateException("Cannot ask cdsend for non-protein coding");
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
	public CodingRNA getCodingRNA(CharSequence chromosomeSequence) {
		return new CodingRNAImpl(chromosomeSequence);
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
