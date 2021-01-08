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

import java.util.stream.IntStream;

import htsjdk.samtools.util.Interval;


public interface TranscriptInterval extends StrandedLocatable {
public Transcript getTranscript();
@Override
default String getContig() {
	return getTranscript().getContig();
}

/** return true if the segment contains the genomic position 1 */
public default boolean contains(final int genomic1) {
	return getStart()<=genomic1 && genomic1 <= getEnd();
}

@Override
public default char getStrand() { return getTranscript().getStrand();}
/** return a name for this interval */
public String getName();

/** return 1-based genomic coordinates from 5' to 3' of the genomic reference */
public default IntStream getGenomicIndexesStream() {
	return IntStream.range(this.getStart(), this.getEnd()+1);
	}

/** convert to an interval */
public default Interval toInterval() {
	return new Interval(
			this.getContig(),
			this.getStart(),
			this.getEnd(),
			this.isNegativeStrand(),
			this.getName()
			);
	}

}
