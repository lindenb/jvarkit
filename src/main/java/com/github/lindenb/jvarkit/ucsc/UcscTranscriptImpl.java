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
	public int hashCode() {
		int i = getTranscriptId().hashCode();
		i= i*31 + getContig().hashCode();
		i= i*31 + Integer.hashCode(txStart);
		i= i*31 + Integer.hashCode(txEnd);
		return i;
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
