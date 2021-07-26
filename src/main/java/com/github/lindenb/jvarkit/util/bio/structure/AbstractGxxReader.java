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
package com.github.lindenb.jvarkit.util.bio.structure;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Locatable;

/** abstract GTF or GFF3 reader */
public abstract class AbstractGxxReader  implements Closeable {
	private static final Logger LOG = Logger.build(AbstractGxxReader.class).make();

	
	
	/** exon coordinate */
	protected static class Coords {
		int start;
		int end;
		}

	
	protected static class TranscriptImpl implements Transcript
		{
		GeneImpl gene;
		String transcript_id;
		int txStart;
		int txEnd;
		StartCodonImpl codon_start= null;
		StopCodonImpl codon_end = null;
		char strand = '?';
		int exonStarts[]=null;
		int exonEnds[]=null;
		final Map<String,String> properties = new HashMap<String, String>();
		boolean coding = false;
		boolean saw_cds_flag = false;
		
		protected abstract class AbstractCodonImpl implements Codon {
			final int pos[]=new int[] {-1,-1,-1};
			
			private void _visit(int loc) {
				int i=0;
				while( i< pos.length) {
					if(pos[i]==-1) break;
					i++;
					}
				if(i>= pos.length) throw new IllegalArgumentException();
				pos[i]=loc;
				}
			
			void visit(Locatable loc) {
				int L=loc.getLengthOnReference();
				if(L>3) throw new IllegalStateException("large codon for "+getName()+" from "+loc);
				for(int i=0;i< L ;i++)
					{
					_visit(loc.getStart()+i);
					}
				Arrays.sort(pos);
				}
			
			private int at(int idx) {
				final int p = pos[idx];
				if(p<1) throw new IllegalStateException("incomplete codon for "+getName()+" "+getTranscript()+ " "+Arrays.toString(this.pos));
				return p;
				}
			
			@Override
			public int getStart() {
				return at(0);
				}
			
			@Override
			public int getMiddle() {
				return at(1);
				}
			
			@Override
			public int getEnd() {
				return at(2);
				}
			
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			}
		
		protected class StartCodonImpl extends AbstractCodonImpl {
			@Override
			public String getName() {
				return "Start Codon";
				}
			}
		
		protected class StopCodonImpl extends AbstractCodonImpl {
			@Override
			public String getName() {
				return "Stop Codon";
				}
			}
		
		protected abstract class AbstractUTR implements UTR {
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int hashCode() {
				return (getTranscript().getId().hashCode() * 31 +this.getStart())*31 + getEnd();
				}

			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof AbstractUTR)) return false;
				final AbstractUTR other = AbstractUTR.class.cast(obj);
				return other.getTranscript().equals(this.getTranscript()) &&
						other.getStart() == this.getStart() &&
						other.getEnd() == this.getEnd()
						;
				}
			}
		
		protected class UTR5Impl extends AbstractUTR {
			
			
			@Override
			public int getStart() {
				return getTranscript().getTxStart();
				}
			public int getEnd() {
				final Optional<Codon> codon;
				if( isPositiveStrand())
					{
					codon= getTranscript().getCodonStart();
					}
				else 
					{
					codon = getTranscript().getCodonStop();
					}
				
				if(!codon.isPresent()) throw new IllegalStateException("no codon for UTR5' for "+getTranscript().getId());
				return codon.get().getStart() - 1 /* one base before the codon */;
				}
			@Override
			public String getName() {
				return (getTranscript().isNegativeStrand()? "3":"5")+"' UTR of "+getTranscript().getId();
				}
			@Override
			public String toString() {
				return getName();
				}
			}

		protected class UTR3Impl extends AbstractUTR {

			@Override
			public int getStart() {
				final Optional<Codon> codon;
				if( isPositiveStrand())
					{
					codon= getTranscript().getCodonStop();
					}
				else 
					{
					codon = getTranscript().getCodonStart();
					}
				if(!codon.isPresent()) throw new IllegalStateException("no codon for UTR3' for "+getTranscript().getId());
				return codon.get().getEnd() + 1 /* one base after the codon */;
				}
			
			public int getEnd() {
				return getTranscript().getTxEnd();
				}
			@Override
			public String getName() {
				return (getTranscript().isNegativeStrand()? "5":"3")+"' UTR of "+getTranscript().getId();
				}
			@Override
			public String toString() {
				return getName();
				}
			}
		protected class CdsImpl implements Cds
			{
			private final int exon_index;
			CdsImpl(int exon_index) {
				this.exon_index = exon_index;
				}
			@Override
			public Exon getExon() {
				return getTranscript().getExon(this.exon_index);
				}
			@Override
			public int getStart() {
				return Math.max(
						getTranscript().getExonStart(this.exon_index),
						getTranscript().getLeftmostCodon().get().getStart()
						);
				}
			@Override
			public int getEnd() {
				return Math.min(
						getTranscript().getExonEnd(this.exon_index),
						getTranscript().getRightmostCodon().get().getEnd()
						);
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.exon_index;
				}
			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof CdsImpl)) return false;
				final CdsImpl c = CdsImpl.class.cast(obj);
				return c.getTranscript().equals(this.getTranscript()) &&
						c.exon_index == this.exon_index
						;
				}
			
			
			@Override
			public String getName() {
				return "CDS."+getExon().getName();
				}
			@Override
			public String toString() {
				return getName();
				}
			}
		
		protected class ExonImpl implements Exon
			{
			private final int index0;
			ExonImpl(final int index0) {
				this.index0 = index0;
				}
			@Override
			public int getIndex() {
				return this.index0;
				}
			@Override
			public int getStart() {
				return getTranscript().getExonStart(this.index0);
				}
			@Override
			public int getEnd() {
				return getTranscript().getExonEnd(this.index0);
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.index0;
				}

			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof ExonImpl)) return false;
				final ExonImpl ex = ExonImpl.class.cast(obj);
				return ex.getTranscript().equals(this.getTranscript()) &&
						ex.getStart() == this.getStart() &&
						ex.getEnd() == this.getEnd()
						;
				}
			
			@Override
			public String getName() {
				return TranscriptImpl.this.transcript_id+ ".Exon"+
							(isNegativeStrand()?
							getExonCount()-this.index0
							:1+this.index0
							);
				}

			
			@Override
			public String toString() {
				return "Exon "+getContig()+":"+getStart()+"-"+getEnd();
				}
			}
		
		protected class IntronImpl implements Intron
			{
			private final int index0;
			IntronImpl(final int index0) {
				this.index0 = index0;
				}
			@Override
			public Transcript getTranscript() {
				return TranscriptImpl.this;
				}
			@Override
			public int getStart() {
				return TranscriptImpl.this.exonEnds[this.index0] + 1;
				}
			@Override
			public int getEnd() {
				return TranscriptImpl.this.exonStarts[this.index0+1] -1;
				}
			
			@Override
			public int hashCode() {
				return TranscriptImpl.this.transcript_id.hashCode() * 31 +this.index0;
				}
			
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof IntronImpl)) return false;
				final IntronImpl ex = IntronImpl.class.cast(obj);
				return ex.getTranscript().equals(this.getTranscript()) &&
						ex.getStart() == this.getStart() &&
						ex.getEnd() == this.getEnd()
						;
				}
			@Override
			public String toString() {
				return "Intron "+getContig()+":"+getStart()+"-"+getEnd();
				}

			@Override
			public String getName() {
				return TranscriptImpl.this.transcript_id+ ".Intron"+
							(isNegativeStrand()?
							getIntronCount()-this.index0
							:1+this.index0
							);
				}
	
			}

		@Override
		public String getId() {
			return this.transcript_id;
			}
		
		@Override
		public Map<String, String> getProperties() {
			return this.properties;
			}
		
		@Override
		public int hashCode() {
			return this.transcript_id.hashCode();
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof TranscriptImpl)) return false;
			final TranscriptImpl tr = TranscriptImpl.class.cast(obj);
			return tr.transcript_id.equals(this.transcript_id);
			}
		
		@Override
		public String getContig() {
			return getGene().getContig();
			}
		@Override
		public int getTxStart() {
			return txStart;
			}
		@Override
		public int getTxEnd() {
			return txEnd;
			}
		
		
		@Override
		public int getStart() {
			return this.getTxStart();
			}
		@Override
		public int getEnd() {
			return this.getTxEnd();
			}
		@Override
		public Gene getGene() {
			return this.gene;
		}
		@Override
		public Optional<Codon> getCodonStart() {
			return Optional.ofNullable(this.codon_start);
		}
		
		@Override
		public boolean hasCodonStartDefined() {
			return this.codon_start!=null;
			}
		
		@Override
		public boolean hasCodonStopDefined() {
			return this.codon_end!=null;
			}
		
		@Override
		public Optional<Codon> getCodonStop() {
			return Optional.ofNullable(this.codon_end);
		}
		
		@Override
		public int getExonStart(int index0) {
			return this.exonStarts[index0];
			}
		@Override
		public int getExonEnd(int index0) {
			return this.exonEnds[index0];
			}
		
		@Override
		public int getExonCount() {
			return this.exonStarts.length;
			}
		
		@Override
		public Exon getExon(final int index0) {
			return new ExonImpl(index0);
			}
		
		@Override
		public List<Exon> getExons() {
			return IntStream.range(0, this.getExonCount()).
					mapToObj(T->getExon(T)).
					collect(Collectors.toList());
			}
		
		@Override
	    public int getIntronCount() {
	    	return this.getExonCount()-1;
	    	}
		@Override
		public List<Intron> getIntrons() {
			return IntStream.range(0, this.getIntronCount()).
					mapToObj(T->getIntron(T)).
					collect(Collectors.toList());
			}
		
		@Override
		public Intron getIntron(int index0) {
			return new IntronImpl(index0);
			}
		@Override
		public Optional<UTR> getUTR5() {
			if(isNonCoding()) return Optional.empty();
			if(isPositiveStrand() && hasCodonStartDefined()) {
				if(this.getStart()>= this.getCodonStart().get().getStart()) return Optional.empty();
				return Optional.of(new UTR5Impl());
				}
			if(isNegativeStrand() && hasCodonStopDefined()) {
				if(this.getStart()>= this.getCodonStop().get().getStart()) return Optional.empty();
				return Optional.of(new UTR5Impl());
				}
			return Optional.empty();
			}
		
		@Override
		public Optional<UTR> getUTR3() {
			if(isNonCoding()) return Optional.empty();
			if(isPositiveStrand() && hasCodonStopDefined()) {
				if(this.getCodonStop().get().getEnd()>= this.getEnd()) return Optional.empty();
				return Optional.of(new UTR3Impl());
				}
			if(isNegativeStrand() && hasCodonStartDefined()) {
				if(this.getCodonStart().get().getEnd()>= this.getEnd()) return Optional.empty();
				return Optional.of(new UTR3Impl());
				}
			return Optional.empty();
			}

		@Override
		public boolean isCoding() {
			return this.coding;
			}
		@Override
		public boolean isNonCoding() {
			return !this.coding;
			}
		
		@Override
		public char getStrand() {
			return this.strand;
			}
		@Override
		public List<Cds> getAllCds() {
			if(!isCoding()) return Collections.emptyList();
			if(!hasCodonStartDefined()) throw new IllegalStateException("transcript has not start codon  defined");
			if(!hasCodonStopDefined()) throw new IllegalStateException("transcript has not end codon defined");
			final List<Cds> L = new ArrayList<>(this.getExonCount());
			for(int i=0;i< this.getExonCount();i++)
				{
				if(this.getExonStart(i) > this.getRightmostCodon().get().getEnd()) break;
				if(this.getExonEnd(i) < this.getLeftmostCodon().get().getStart()) continue;
				L.add(new CdsImpl(i));
				}
			return L;
			}
		
		@Override
		public String toString() {
			return this.transcript_id+" "+getContig()+":"+getStart()+"-"+getEnd();
			}
		}
	
	protected static class GeneImpl implements Gene
		{
		String gene_id;
		String contig;
		int start;
		int end;
		char strand;
		final List<Transcript> transcripts= new ArrayList<>();
		final Map<String,String> properties = new HashMap<String, String>();
		
		
		@Override
		public String getId() {
			return this.gene_id;
			}
		
		@Override
		public List<Transcript> getTranscripts() {
			return transcripts;
			}
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		@Override
		public Map<String, String> getProperties() {
			return properties;
			}
		@Override
		public char getStrand() {
			return this.strand;
			}
		@Override
		public int hashCode() {
			return this.gene_id.hashCode();
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof GeneImpl)) return false;
			final GeneImpl tr = GeneImpl.class.cast(obj);
			return tr.gene_id.equals(this.gene_id);
			}

		@Override
		public String toString() {
			return this.gene_id+" "+getContig()+":"+getStart()+"-"+getEnd();
			}
		
		}
}
