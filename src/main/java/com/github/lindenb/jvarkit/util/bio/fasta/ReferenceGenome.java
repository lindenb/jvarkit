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

import java.io.Closeable;
import java.util.function.Function;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public interface ReferenceGenome extends
	Closeable,
	Function<String, ReferenceContig> {
public SAMSequenceDictionary getDictionary();
public ReferenceContig getContig(final String contigName);

/** convert contig name using aliases of the SAMSequenceDictionary, return null if there is no alias */
public default String convertContig(final String srcName)
	{
	final SAMSequenceRecord rec = getDictionary().getSequence(srcName);
	return rec==null?null:rec.getSequenceName();
	}

/** returns a ReferenceContig by tid or null if it doesn't exist */
public default ReferenceContig getContig(final int id) {
	return getContig(getDictionary().getSequence(id).getSequenceName());
}

/** return true if there is a contig with this name or alias */
public default boolean hasContig(final String name)
	{
	return getDictionary().getSequence(name)!=null;
	}

/** return the dictionary size */
public default int size() { return getDictionary().size();}

/** return true if there is the dictionary is empty */
public default boolean isEmpty() { return getDictionary().isEmpty();}

public String getSource();

@Override
public default ReferenceContig apply(final String contigName) {
		return getContig(contigName);
	}
}
