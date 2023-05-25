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

import java.util.AbstractList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;

public class UcscTranscript implements Locatable {
	private String contig;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	private char strand;
	private int[] exonStarts;
	private int[] exonEnds;
	private String name;
	
	public UcscTranscript(final String s) {
		this(CharSplitter.TAB.split(s));
		}
	
	public UcscTranscript(final String[] tokens) {

		}
	private abstract class Coordinates {
		private final int[] genomicIndexes0;
		Coordinates(List<Locatable> L) {
			genomicIndexes0 = new int[L.stream().mapToInt(X->X.getLengthOnReference()).sum()];
			int i=0;
			Collections.sort(L,(A,B)->Integer.compare(A.getStart(), B.getStart()));
			for(Locatable loc: L) {
				for(int x=loc.getStart();x<=loc.getEnd();++x) {
					genomicIndexes0[i++] = x-1;
					}
				}
			}
		public int convertToGenomicIndex0(int mRna0) {
			if(isPositiveStrand()) {
				return genomicIndexes0[mRna0];
				}
			else
				{
				return genomicIndexes0[(genomicIndexes0.length-1) - mRna0];
				}
			}
		public int size() {
			return genomicIndexes0.length;
			}
		}
	/** ===================================================================*/
	public abstract class Component implements Locatable {
		public UcscTranscript getTranscript() {
			return UcscTranscript.this;
			}
		
		@Override
		public int hashCode() {
			int i = getContig().hashCode();
			i= i*31 + getTranscript().getName().hashCode();
			i= i*31 + Integer.hashCode(getStart());
			i= i*31 + Integer.hashCode(getEnd());
			return i;
			}
		
		public abstract String getName();
		
		public Strand getStrand() {
			return UcscTranscript.this.getStrand();
			}
		
		public boolean isPositiveStrand() {
			return UcscTranscript.this.isPositiveStrand();
		}
		
		public boolean isNegativeStrand() {
			return UcscTranscript.this.isNegativeStrand();
		}
		
		@Override
		public String getContig() {
			return getTranscript().getContig();
			}
		@Override
		public String toString() {
			return getClass().getSimpleName()+":"+getTranscript().getName()+":"+getContig()+":"+getStart()+"-"+getEnd();
			}	
		}
	/** ===================================================================*/
	public class Exon extends Component {
		private int exon_index;
		Exon(int exon_index) {
			this.exon_index = exon_index;
			}
		
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof Exon)) return false;
			final Exon o = Exon.class.cast(obj);
			return exon_index==o.exon_index &&
					contigsMatch(o) && 
					getTranscript().getName().equals(o.getTranscript().getName());
			}
		
		@Override
		public int getStart() {
			return UcscTranscript.this.exonStarts[this.exon_index] +1 ;
			}
		@Override
		public int getEnd() {
			return UcscTranscript.this.exonEnds[this.exon_index] ;
			}
		@Override
		public String getName() {
			return getTranscript().getName()+".Exon"+(isPositiveStrand()?exon_index+1:getExonCount()-exon_index);
			}
		}
	/** ===================================================================*/
	public class Intron extends Component {
		private int intron_index;
		private Intron(int intron_index) {
			this.intron_index = intron_index;
			}
		@Override
		public int getStart() {
			return UcscTranscript.this.exonEnds[this.intron_index]  ;
			}
		@Override
		public int getEnd() {
			return UcscTranscript.this.exonStarts[this.intron_index+1] ;
			}
		@Override
		public String getName() {
			return getTranscript().getName()+".Intron"+(isPositiveStrand()?intron_index+1:getIntronCount()-intron_index);
			}
		
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof Intron)) return false;
			final Intron o = Intron.class.cast(obj);
			return intron_index==o.intron_index && 
					contigsMatch(o) && 
					getTranscript().getName().equals(o.getTranscript().getName());
			}
		}
	/** ===================================================================*/
	public abstract class ExonComponent extends Component {
		private final Exon exon;
		private ExonComponent(final Exon exon) {
			this.exon = exon;
			}
		public final Exon getExon() {
			return this.exon;
			}
		@Override
		public String getName() {
			return getTranscript().getName()+"."+getClass().getSimpleName() + "."+getStart()+"-"+getEnd();
			}
		}
	
	/** ===================================================================*/

	public class CDS extends ExonComponent {
		private CDS(final Exon exon) {
			super(exon);
			}
		@Override
		public int getStart() {
			return Math.max(getExon().getStart(), cdsStart+1);
			}
		@Override
		public int getEnd() {
			return Math.min(getExon().getEnd(), cdsEnd);
			}
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof CDS)) return false;
			final CDS o = CDS.class.cast(obj);
			return getStart()==o.getStart() && 
					getEnd()==o.getEnd() && 
					contigsMatch(o) && 
					getTranscript().getName().equals(o.getTranscript().getName());
			}
		}
	/** ===================================================================*/
	public abstract class UTR extends ExonComponent {
		protected UTR(final Exon exon) {
			super(exon);
			}
		
		}
	/** ===================================================================*/
	public class UTR5 extends UTR {
		private UTR5(final Exon exon) {
			super(exon);
			}
		@Override
		public int getStart() {
			if(isPositiveStrand()) {
				return getExon().getStart();
				}
			else
				{
				return Math.max(getExon().getEnd(), cdsEnd);
				}
			}
		@Override
		public int getEnd() {
			if(isPositiveStrand()) {
				return Math.min(getExon().getEnd(), cdsStart /* no +1, base before */);
				}
			else
				{
				return getExon().getEnd();
				}
			}
		}
	/** ===================================================================*/
	public class UTR3 extends UTR {
		private UTR3(final Exon exon) {
			super(exon);
			}
		@Override
		public int getStart() {
			if(isPositiveStrand()) {
				return Math.max(getExon().getStart(), cdsEnd);
				}
			else
				{
				return getExon().getStart();
				}
			}
		@Override
		public int getEnd() {
			if(isPositiveStrand()) {
				return getExon().getEnd();
				}
			else
				{
				return Math.min(getExon().getEnd(), cdsStart /* no +1, base before */);
				}
			}

		}
	/** ===================================================================*/
	
	
	
	public Strand getStrand() {
		return Strand.decode(this.strand);
		}
	
	public boolean isPositiveStrand() {
		return this.strand=='+';
	}
	
	public boolean isNegativeStrand() {
		return this.strand=='-';
	}
	
	@Override
	public String getContig() {
		return this.contig;
		}
	
	@Override
	public int getStart() {
		return getTxStart();
		}
	
	@Override
	public int getEnd() {
		return getTxEnd();
		}
	
	public boolean isProteinCoding() {
		return this.cdsStart < this.cdsEnd;
		}
	
	public OptionalInt getCdsStart() {
		return isProteinCoding()?OptionalInt.of(this.cdsStart+1):OptionalInt.empty();
		}
	
	public OptionalInt getCdsEnd() {
		return isProteinCoding()?OptionalInt.of(this.cdsEnd):OptionalInt.empty();
		}

	
	
	public String getName() {
		return name;
		}
	
	@Override
	public int hashCode() {
		int i = getName().hashCode();
		i= i*31 + getContig().hashCode();
		i= i*31 + Integer.hashCode(txStart);
		i= i*31 + Integer.hashCode(txEnd);
		return i;
		}
	
	public int getTxStart() {
		return txStart+1;
		}
	public int getTxEnd() {
		return txEnd;
		}
	
	
	public int getIntronCount() {
		return Math.max(0, getExonCount()-1);
		}
	
	public int getExonCount() {
		return this.exonStarts.length;
		}
	
	public Intron getIntron(int idx) {
		return new Intron(idx);
		}
	
	public Exon getExon(int idx) {
		return new Exon(idx);
		}
	
	public List<Exon> getExons() {
		return new AbstractList<Exon>() {
			@Override
			public Exon get(int idx) {
				return getExon(idx);
				}
			@Override
			public int size() {
				return getExonCount();
				}
			};
		}
	public List<Intron> getIntrons() {
		if(getIntronCount()==0) return Collections.emptyList();
		return new AbstractList<Intron>() {
			@Override
			public Intron get(int idx) {
				return getIntron(idx);
				}
			@Override
			public int size() {
				return getIntronCount();
				}
			};
		}
	
	public List<CDS> getCDSs() {
		if(!isProteinCoding()) return Collections.emptyList();
		return getExons().stream().
				filter(EX->CoordMath.overlaps(EX.getStart(), EX.getEnd(), this.cdsStart+1, this.cdsEnd)).
				map(EX->new CDS(EX)).
				collect(Collectors.toList());
		}

	
	@Override
	public String toString() {
		StringBuilder sb= new StringBuilder();
		return sb.toString();
		}
	}
