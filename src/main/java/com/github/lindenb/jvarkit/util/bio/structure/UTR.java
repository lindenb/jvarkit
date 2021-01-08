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

import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

public interface UTR extends TranscriptInterval {
	
/** return the intervals for this utr: the list may have more than one interval if the UTR is spliced */
public default List<Locatable> getIntervals() {
	return getTranscript().
		getExons().
		stream().
		filter(E->E.overlaps(this)).
		map(E->new Interval(
				getContig(),
				Math.max(E.getStart(), this.getStart()),
				Math.min(E.getEnd(), this.getEnd()),
				isNegativeStrand(),
				"UTR in "+E.getName()
				)
				).
		collect(Collectors.toList());
		}

/** return true wether this urt is spliced between two exons */
public default boolean isSpliced() {
	return getIntervals().size()>1;
	}
}
