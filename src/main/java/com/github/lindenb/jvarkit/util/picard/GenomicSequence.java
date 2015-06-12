/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import htsjdk.samtools.SAMSequenceRecord;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;

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
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private SAMSequenceRecord samSequenceRecord;
	private byte buffer[]=null;
	private int buffer_pos=-1;
	private int half_buffer_capacity=1000000;
	public GenomicSequence(IndexedFastaSequenceFile indexedFastaSequenceFile ,String chrom)
		{	
		this.indexedFastaSequenceFile=indexedFastaSequenceFile;
		if(this.indexedFastaSequenceFile==null) throw new NullPointerException("IndexedFastaSequenceFile is null");
		if(this.indexedFastaSequenceFile.getSequenceDictionary()==null)
			{
			throw new NullPointerException("no sequence dictionary in the reference. Use picard CreateSequenceDictionary to index the sequence https://broadinstitute.github.io/picard/command-line-overview.html.");
			}
		this.samSequenceRecord=this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(chrom);
		if(this.samSequenceRecord==null) throw new IllegalArgumentException("not chromosome "+chrom+" in reference.");
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
	}
