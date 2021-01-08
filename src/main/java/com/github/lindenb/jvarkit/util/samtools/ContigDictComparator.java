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
package com.github.lindenb.jvarkit.util.samtools;

import java.util.Comparator;

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

/**
 * String comparator using a SAMSequenceDictionary
 */
public class ContigDictComparator implements Comparator<String> {

private final SAMSequenceDictionary dict;

public ContigDictComparator(final SAMSequenceDictionary dict) {
	this.dict = dict;
	}
protected int convertToTid(final String name) {
	final int tid = this.dict.getSequenceIndex(name);
	if( tid == -1) throw new JvarkitException.ContigNotFoundInDictionary(name, this.dict);
	return tid;
	}
@Override
public int compare(final String o1,final String o2) {
	return Integer.compare( convertToTid(o1) , convertToTid(o2));
	}

/** creates a generic comparator for Locatable using this ContigDictComparator to compare the contigs */
public <T extends Locatable> Comparator<T> createLocatableComparator() {
	return (A,B) ->{
		int i = ContigDictComparator.this.compare(A.getContig(),B.getContig());
		if(i!=0) return i;
		i =Integer.compare(A.getStart(), B.getStart());
		if(i!=0) return i;
		return Integer.compare(A.getEnd(), B.getEnd());
		};
	}

}
