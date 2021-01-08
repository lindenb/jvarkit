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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.function.Function;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;

import htsjdk.samtools.util.Locatable;

public class RNASequenceFactory {
	
	@FunctionalInterface
	public static interface BaseConverter {
		public char transform(char b);
	}
	
	private Function<String,CharSequence> contig2genomicSequence = S->{throw new IllegalArgumentException();};
	private BaseConverter baseConverter = B->Character.toUpperCase(B);
	
	public RNASequenceFactory() {		
	}
	
	public RNASequenceFactory setContigToGenomicSequence(final Function<String, CharSequence> contig2genomicSequence) {
		this.contig2genomicSequence = contig2genomicSequence;
		return this;
		}
	
	public RNASequenceFactory setBaseConverter(final BaseConverter baseConverter) {
		this.baseConverter = baseConverter;
		return this;
		}
	
	private CharSequence getRequiredGenomicSequence(final Transcript tr) {
		final CharSequence chrom = this.contig2genomicSequence.apply(tr.getContig());
		if(chrom==null) throw new IllegalArgumentException("no contig available for "+tr.getContig());
		return chrom;
	}
	
	public RNASequence getMessengerRNA(final Transcript tr) {
		return new BasicRNASequence(tr, getRequiredGenomicSequence(tr),tr.getExons());
	}
	
	public RNASequence getCodingRNA(final Transcript tr) {
		return new BasicRNASequence(tr, getRequiredGenomicSequence(tr),tr.getAllCds());
	}
	
	private class BasicRNASequence extends AbstractCharSequence implements RNASequence {
		final Transcript transcript;
		final int length;
		final int rna02genomic0[];
		final byte bases[];
		final Map<Integer, Integer> genomic0torna0;
		protected BasicRNASequence(final Transcript transcript,final CharSequence contigSeq,
				final List<? extends Locatable> intervals) {
			this.transcript = transcript;			
			this.length = intervals.stream().mapToInt(R->R.getLengthOnReference()).sum();
			this.rna02genomic0 = new int[this.length];
			this.bases = new byte[this.length];
			this.genomic0torna0 = new HashMap<>(this.length);
			
			if(transcript.isPositiveStrand()) {
				int x=0;
				for(Locatable r:intervals) {
					for(int i=r.getStart();i<=r.getEnd();++i) {
						final int g0 = i-1;
						this.rna02genomic0[x] = g0;
						this.genomic0torna0.put(g0,x);
						this.bases[x]  = (byte)baseConverter.transform(contigSeq.charAt(g0));
						++x;
						}
					}
				}
			else
				{
				int x=0;
				for(int r1 = intervals.size()-1;r1>=0;r1--) {
					final Locatable r = intervals.get(r1);
					for(int i=r.getEnd();i>=r.getStart();--i) {
						final int g0 = i-1;
						this.rna02genomic0[x] = g0;
						this.genomic0torna0.put(g0,x);
						this.bases[x]  = (byte)AcidNucleics.complement(baseConverter.transform(contigSeq.charAt(g0)));
						++x;
						}
					}
				}
			}
		
		@Override
		public char charAt(int index) {
			return (char)this.bases[index];
			}
		
		@Override
		public OptionalInt convertGenomic0ToRnaIndex0(final int g0) {
			final Integer p0 = this.genomic0torna0.get(g0);
			return p0==null?OptionalInt.empty():OptionalInt.of(p0.intValue());
			}
		
		@Override
		public int convertRnaIndex0ToGenomic0(int rna0) {
			return this.rna02genomic0[rna0];
			}
		
		@Override
		public int length() {
			return this.length;
			}
		@Override
		public Transcript getTranscript() {
			return this.transcript;
			}
		
	}
	
}
