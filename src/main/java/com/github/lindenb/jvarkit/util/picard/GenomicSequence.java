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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.picard;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import java.util.OptionalInt;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

/**
 * 
 * implementation of java.lang.CharSequence for a given
 * chromosome of a picard IndexedFastaSequenceFile
 *
 */
public class GenomicSequence
	extends AbstractCharSequence
	implements ChromosomeSequence
	{
	private final IndexedFastaSequenceFile indexedFastaSequenceFile;
	private final SAMSequenceRecord samSequenceRecord;
	private byte buffer[]=null;
	private int buffer_pos=-1;
	private int half_buffer_capacity=1000000;
	
	public static interface GCPercent extends Locatable
		{
		public int getAllCount();
		public int getGCCount();
		public int getATCount();
		/** return true if getAllCount==0 */
		public boolean isEmpty();
		/** return GC% as double between 0 and 1 . return -1 if interval isEmpty */
		public double getGCPercent();
		/** return GC% as int between 0 and 100 . return -1 if interval isEmpty */
		public int getGCPercentAsInteger();
		/** return GC% as optional */
		public default OptionalInt getOptGCPercent() {
			return isEmpty()?
					OptionalInt.empty():
					OptionalInt.of(getGCPercentAsInteger())
					;
			}

		}
	
	private static class GCPercentImpl
		implements GCPercent
		{
		final String contig;
		final int start1;
		final int end1;
		int count=0;
		int count_gc=0;
		int count_at=0;

		GCPercentImpl(String contig,int s1,int e1) {
			this.contig = contig;
			this.start1=s1;
			this.end1=e1;
			}
		@Override public int getAllCount() { return this.count;}
		@Override public int getGCCount() { return this.count_gc;}
		@Override public int getATCount(){ return this.count_at;}
		@Override
		public boolean isEmpty() { return this.count == 0; }
		@Override
		public double getGCPercent() {
			return ( this.count==0? -1.0:(this.count_gc/(double)this.count));
			}
		@Override
		public int getGCPercentAsInteger() {
			return ( this.count==0? -1:(int)(getGCPercent()*100.0));
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + count;
			result = prime * result + count_at;
			result = prime * result + count_gc;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null || !(obj instanceof GCPercentImpl)) {
				return false;
				}
			final GCPercentImpl other = (GCPercentImpl) obj;
			return	this.count==other.count &&
					this.count_at==other.count_at &&
					this.count_gc==other.count_gc;
			}
		@Override
		public String toString() {
			return "gc_percent("+this.count_gc+"/"+this.count+")="+this.getGCPercent();
			}
		@Override
		public String getContig() {
			return this.contig;
		}
		@Override
		public int getStart() {
			return this.start1;
		}
		@Override
		public int getEnd() {
			return this.end1;
		}
		}
	
	public GenomicSequence(final IndexedFastaSequenceFile indexedFastaSequenceFile ,final String chrom)
		{	
		this.indexedFastaSequenceFile=indexedFastaSequenceFile;
		if(this.indexedFastaSequenceFile==null) throw new NullPointerException("IndexedFastaSequenceFile is null");
		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
		this.samSequenceRecord= dict.getSequence(chrom);
		if(this.samSequenceRecord==null) throw new JvarkitException.ContigNotFoundInDictionary(chrom,dict);
		}
	
	public SAMSequenceRecord getSAMSequenceRecord()
		{
		return samSequenceRecord;
		}
	
	@Override
	public int hashCode()
		{
		return getSAMSequenceRecord().hashCode();
		}
	
	/** get the chromosome name of that genomic sequence */
	@Override
	public String getChrom()
		{
		return getSAMSequenceRecord().getSequenceName();
		}
	
	@Override
	public int length()
		{
		return getSAMSequenceRecord().getSequenceLength();
		}
	
	@Override
	public char charAt(int index0)
		{
		if(index0 >= length())
			{
			throw new IndexOutOfBoundsException("index:"+index0);
			}
		if(buffer!=null && index0>=buffer_pos && index0-buffer_pos < buffer.length)
			{
			return (char)buffer[index0-buffer_pos];
			}
		int minStart=Math.max(0, index0-half_buffer_capacity);
		int maxEnd=Math.min(minStart+2*half_buffer_capacity,this.length());
		this.buffer=this.indexedFastaSequenceFile.getSubsequenceAt(
				getChrom(),
				minStart+1,
				maxEnd).getBases();
		this.buffer_pos=minStart;
		return (char)buffer[index0-minStart];
		}
	
	/** return GC% between start (inclusive, 0 based) and end (exclusive)) */
	public GCPercent getGCPercent(int start,int end) {
		final int L=this.length();
		final GCPercentImpl gcp = new GCPercentImpl(
				this.getChrom(),
				start+1,
				Math.min(end, L)
				);
		for(int i=start;i< end && i< L;++i) {
			gcp.count++;
			switch(this.charAt(i)) {
				case 'c': case 'C':
				case 'g': case 'G':
				case 's': case 'S':gcp.count_gc++; break;
				case 'a': case 'A':
				case 't': case 'T':
				case 'w': case 'W':gcp.count_at++; break;
				}
			}
		return gcp;
		}
	
	
	}
