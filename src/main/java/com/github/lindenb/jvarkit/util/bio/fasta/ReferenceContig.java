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

*/
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.util.OptionalDouble;
import java.util.OptionalInt;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;

public interface ReferenceContig
	extends CharSequence,Locatable{

	
	public static interface GCPercent extends Locatable
		{
		public int getAllCount();
		public int getGCCount();
		public int getATCount();
		/** return true if getAllCount==0 */
		public boolean isEmpty();
		/** return GC% as double between 0 and 1  */
		public OptionalDouble getGCPercent();
		/** return GC% as int between 0 and 100 . return -1 if interval isEmpty */
		public OptionalInt getGCPercentAsInteger();
		}
	
	
public SAMSequenceRecord getSAMSequenceRecord();

/** return true if contig is compatible with 'name' */
public default boolean hasName(final String name) {
	return getContig().equals(name);
	}

@Override
public	default String getContig() {
		return getSAMSequenceRecord().getSequenceName();
		}

@Override
public default int length() {
	return getSAMSequenceRecord().getSequenceLength();
	}

@Override 
public default int getStart() {
	return 1;
	}
@Override 
public default int getEnd() {
	return this.length();
	}

/** return GC%
 * 
 * @param start 0 based start pos
 * @param end 0 based end pos
 * @return a GC% object
 */
public GCPercent getGCPercent(int start0,int end0);

}
