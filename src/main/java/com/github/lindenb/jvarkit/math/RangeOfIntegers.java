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
package com.github.lindenb.jvarkit.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.lang.CharSplitter;

public class RangeOfIntegers {
public static final String OPT_DESC="A 'range of integers' is a list of integers in ascending order separated with semicolons. ";
/** generate RangeOfIntegers from a String with values separated with semicolons */
public static class StringConverter 
	implements com.beust.jcommander.IStringConverter<RangeOfIntegers>
	{
	@Override
	public RangeOfIntegers convert(final String s) {
		return new RangeOfIntegers(s);
		}
	}
	
public static interface Range
	extends Comparable<Range>
	{
	public Integer getMinInclusive();
	public Integer getMaxExclusive();
	public boolean contains(int value);
	}

private static class RangeImpl implements Range
	{
	final Integer minIncl;
	final Integer maxExcl;
	RangeImpl(final Integer minIncl,final Integer maxExcl) {
		this.minIncl=minIncl;
		this.maxExcl=maxExcl;
		}
	@Override
	public Integer getMinInclusive() {
		return minIncl;
		}
	@Override
	public Integer getMaxExclusive() {
		return maxExcl;
		}
	
	@Override
	public int compareTo(final Range o) {
		if(this.getMinInclusive()==null) return -1;
		if(o.getMinInclusive()==null) return 1;
		return getMinInclusive().compareTo(o.getMinInclusive());
		}
	@Override
	public boolean contains(final int value) {
		if(this.minIncl!=null && value < this.minIncl) return false;
		if(this.maxExcl!=null && value >= this.maxExcl) return false;
		return true;
		}
	@Override
	public String toString() {
		if(minIncl!=null && maxExcl!=null && this.minIncl+1 == this.maxExcl)
			{
			return String.valueOf(this.minIncl);
			}
		
		return "[" + 
				(minIncl==null?"-Inf":String.valueOf(this.minIncl))+ "/"+ 
				(maxExcl==null?"Inf":String.valueOf(this.maxExcl)) +
				"[";
		}
	}


private final List<Range> ranges;

/** generate RangeOfIntegers from a String with values separated with semicolons */
public RangeOfIntegers(final String s) {
	this(Arrays.asList(CharSplitter.SEMICOLON.split(s)).stream().
			filter(S->!S.trim().isEmpty()).
			mapToInt(S->Integer.valueOf(S)).
			toArray());
	}

public RangeOfIntegers(final int...array) {
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

public Range getRange(final int value) {
	for(final Range r:getRanges())
		{
		if(r.contains(value)) return r;
		}
	throw new IllegalStateException("cannot get range ??" +value);
	}
@Override
public String toString() {
		return getRanges().toString();
	}
}
