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
package com.github.lindenb.jvarkit.fastq;

import java.util.AbstractList;
import java.util.List;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;

/** a pair of FastqRecord */
public interface FastqRecordPair extends List<FastqRecord> {
/** get first of pair */
public FastqRecord getFirstInPair();
/** get second of pair */
public FastqRecord getSecondInPair();
/** check names R1 and R2 are valid */
public boolean isValidName();
/** assert names R1 and R2 are valid  and return this*/
public default FastqRecordPair validateName() {
	if(!isValidName()) throw new IllegalArgumentException("Bad Fastq R1 and R2 names: "
			+ getFirstInPair().getReadName()+" "
			+ getSecondInPair().getReadName()
			);
	return this;
	}


public static FastqRecordPair of(final FastqRecord fq1, final FastqRecord fq2) {
	return new FastqRecordPairImpl(fq1,fq2);
	}

static class FastqRecordPairImpl extends AbstractList<FastqRecord>  implements FastqRecordPair {
	private final FastqRecord r1;
	private final FastqRecord r2;
	FastqRecordPairImpl(final FastqRecord r1,final FastqRecord r2) {
		this.r1 = r1;
		this.r2 = r2;
		}
	
	@Override
	public FastqRecord getFirstInPair() {
		return this.r1;
		}
	@Override
	public FastqRecord getSecondInPair() {
		return this.r2;
		}
	
	@Override
	public final int size() {
		return 2;
		}
	
	@Override
	public FastqRecord get(final int index) {
		switch(index) {
			case 0: return this.r1;
			case 1: return this.r2;
			default: throw new IndexOutOfBoundsException("0<="+index+"<2");
			}
		}
	
	/** normalize read name */
	private String normalize(String s, int side) {
		int w = s.indexOf(' ');
		if(w==-1) w=s.indexOf('\t');
		if(w!=-1) s=s.substring(0,w);
		if(
			(side==0 && s.endsWith(FastqConstants.FIRST_OF_PAIR)) || 
			(side==1 && s.endsWith(FastqConstants.SECOND_OF_PAIR))) {
			s=s.substring(0,s.length()-2);
			}
		return s;
		}
	
	@Override
	public boolean isValidName() {
		String n1 = this.getFirstInPair().getReadName();
		String n2 = this.getSecondInPair().getReadName();
		if(n1.equals(n2)) return true;
		n1 = normalize(n1,0);
		n2 = normalize(n2,0);
		return n1.equals(n2);
		}
	
	@Override
	public int hashCode() {
		return this.r1.hashCode()*31 + this.r2.hashCode();
		}
	@Override
	public boolean equals(final Object o) {
		if(o==this) return true;
		if(o==null || !(o instanceof FastqRecordPair)) return false;
		final FastqRecordPair other = FastqRecordPair.class.cast(o);
		return this.r1.equals(other.getFirstInPair()) &&
				this.r2.equals(other.getSecondInPair())
				;
		}
	@Override
	public String toString() {
		return "FastqRecordPair("+this.r1+","+this.r2+")";
		}
	}
}
