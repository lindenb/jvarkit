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
package com.github.lindenb.jvarkit.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

public class RangeOfDoubles {
	public static final String OPT_DESC="A 'range of double' is a list of floating number in ascending order separated with semicolons. ";

	/** generate RangeOfIntegers from a String with values separated with semicolons */
public static class StringConverter 
	implements com.beust.jcommander.IStringConverter<RangeOfDoubles>
	{
	@Override
	public RangeOfDoubles convert(final String s) {
		return new RangeOfDoubles(s);
		}
	}
	
public static interface Range
	extends Comparable<Range>
	{
	public Double getMinInclusive();
	public Double getMaxExclusive();
	/** return true if getMinInclusive()!=null */
	public default boolean isMinDefined() {
		return getMinInclusive()!=null;
		}
	/** return true if getMaxInclusive()!=null */
	public default boolean isMaxDefined() {
		return getMaxExclusive()!=null;
		}
	public boolean contains(double value);
	}

private static class RangeImpl implements Range
	{
	final Double minIncl;
	final Double maxExcl;
	
	RangeImpl(final Double minIncl,final Double maxExcl) {
		this.minIncl=minIncl;
		this.maxExcl=maxExcl;
		}
	@Override
	public Double getMinInclusive() {
		return minIncl;
		}
	@Override
	public Double getMaxExclusive() {
		return maxExcl;
		}
	
	@Override
	public int hashCode() {
		return (minIncl==null?0:minIncl.hashCode())*31  + (maxExcl==null?0:maxExcl.hashCode());
		}
	public boolean equals(final Object o) {
		if(o==this) return true;
		if(o==null || !(o instanceof RangeImpl)) return false;
		return compareTo(Range.class.cast(o))==0;
		}

	@Override
	public int compareTo(final Range o) {
		if(o==this) return 0;
		if(this.getMinInclusive()==null) return -1;
		if(o.getMinInclusive()==null) return 1;
		return getMinInclusive().compareTo(o.getMinInclusive());
		}
	@Override
	public boolean contains(final double value) {
		if(this.minIncl!=null && value < this.minIncl) return false;
		if(this.maxExcl!=null && value >= this.maxExcl) return false;
		return true;
		}
	@Override
	public String toString() {
		return "[" + 
				(minIncl==null?"-Inf":String.valueOf(this.minIncl))+ " / "+ 
				(maxExcl==null?"Inf":String.valueOf(this.maxExcl)) +
				"[";
		}
	}


private final List<Range> ranges;

/** generate RangeOfIntegers from a String with values separated with semicolons */
public RangeOfDoubles(final String s) {
	this(Arrays.asList(s.split(";")).stream().filter(S->!S.trim().isEmpty()).mapToDouble(S->Double.valueOf(S)).toArray());
	}

public RangeOfDoubles(final double array[]) {
	if(array==null) throw new IllegalArgumentException("array cannot be null");
	if(array.length==0) throw new IllegalArgumentException("array cannot be empty");
	final List<Range> ranges = new ArrayList<>(array.length+2);
	ranges.add(new RangeImpl(null, array[0]));
	for(int i=0;i+1< array.length;++i)
		{
		if(array[i]>=array[i+1])
			{
			throw new IllegalArgumentException("limit should be specified in ascending order :"+Arrays.toString(array));
			}
		ranges.add(new RangeImpl( array[i], array[i+1]));
		}
	ranges.add(new RangeImpl(array[array.length-1],null));
	this.ranges = Collections.unmodifiableList(ranges);
	}

public List<Range> getRanges() {
	return this.ranges;
	}

/** alias of getRanges().stream() */
public Stream<Range> stream() {
	return getRanges().stream();
	}


public Range getRange(final double value) {
	for(final Range r:getRanges())
		{
		if(r.contains(value)) return r;
		}
	throw new IllegalStateException("cannot get range ??" +value);
	}

/**
 * generate a new RangeOfDoubles from a sequence of numbers
 * @param begin start value inclusive
 * @param end end value inclusive
 * @param step shift value
 * @return new RangeOfDouble
 */
public static RangeOfDoubles fromTo(double begin,double end,double step) {
	if(begin>end) throw new IllegalArgumentException(""+begin+">"+end);
	if(step<=0) throw new IllegalArgumentException(""+step+"<=0");
	final List<Double> L = new ArrayList<>();
	while(begin<=end) {
		L.add(begin);
		begin+=step;
		}
	return new RangeOfDoubles(L.stream().mapToDouble(D->D.doubleValue()).toArray());
	}
}
