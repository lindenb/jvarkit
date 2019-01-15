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

import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.iterator.MergingIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.CloseableIterator;

/** a fast ? version of merging samrecorditerator, no merging of read groups
 * assuming sam sequence dictionaries are the same */
@Deprecated // use MergingIterator
public class MergingSamRecordIterator 
	extends MergingIterator<SAMRecord>
	{
	public MergingSamRecordIterator(
			final Comparator<SAMRecord> comparator,
			List<CloseableIterator<SAMRecord>> iterators
			)
		{
		super(comparator,iterators);
		}
	public MergingSamRecordIterator(final List<CloseableIterator<SAMRecord>> iterators)
		{
		this(new SAMRecordCoordinateComparator(),iterators);
		}
	}
