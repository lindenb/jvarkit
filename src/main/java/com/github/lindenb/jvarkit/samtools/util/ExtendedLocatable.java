/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.samtools.util;


import com.github.lindenb.jvarkit.bed.BedInterval;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

public interface ExtendedLocatable extends Locatable, BedInterval {


/** return true if the interval contains this 1-based position */
public default boolean contains(int g1) {
	return this.getStart() <=g1 && g1 <=this.getEnd();
	}

@Override
public default int getBedStart() {
	return getStart() - 1;
	}
@Override
public default int getBedEnd() {
	return getEnd();
	}

/** return contig(tab)(start-1)(tab)(end) */
public default String toBed3() {
	return LocatableUtils.toBed3(this);
	}

/** use StringUtils.niceInt to format start and end  */
public default String toNiceString() {
	return LocatableUtils.toNiceString(this);
	}

@Override
public default Locatable toLocatable() {
	return this;
	}

/** convert this to Interval, useful to insert in IntervalTreeMap */
public default Interval toInterval() {
	if(this instanceof Interval) return Interval.class.cast(this);
	return new Interval(this);
	}
/** if this and other are not on the same contig, throws an exception
 *  if they overlap return 0
 *  otherwise return distance from end to start
 **/
public default int getDistanceTo(final Locatable other) {
	if(!contigsMatch(other)) throw new IllegalArgumentException("not the same contigs");
	if(CoordMath.overlaps(getStart(), getEnd(), other.getStart(), other.getEnd())) return 0;
	if(this.getEnd() < other.getStart()) {
		return CoordMath.getLength(getEnd(), other.getStart());
		}
	else if(other.getEnd() < this.getStart()) {
		return CoordMath.getLength(other.getEnd(), this.getStart());
		}
	else
		{
		throw new IllegalStateException();
		}
	}
}
