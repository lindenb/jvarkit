/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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

import java.util.Comparator;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

public class LocatableUtils extends CoordMath {
	public static final Comparator<Locatable> DEFAULT_COMPARATOR = (A,B)->LocatableUtils.compareTo(A,B);
	
	/** generic comparator on chrom/start/end */
	public static int compareTo(final Locatable A, Locatable B) {
		int i= A.getContig().compareTo(B.getContig());
		if(i!=0) return i;
		i = Integer.compare(A.getStart(),B.getStart());
		if(i!=0) return i;
		i = Integer.compare(A.getEnd(),B.getEnd());
		return i;
		}
	
	
	/** use chrom:start-end  */
	public static String toString(final Locatable loc) {
		return loc.getContig()+":"+ 
				String.valueOf(loc.getStart()) + "-"+ 
				String.valueOf(loc.getEnd())
				;
		}
	
	/** use StringUtils.niceInt to format start and end  */
	public static String toNiceString(final Locatable loc) {
		return loc.getContig()+":"+ 
				StringUtils.niceInt(loc.getStart()) + "-"+ 
				StringUtils.niceInt(loc.getEnd())
				;
		}
	/** use StringUtils.niceInt to format start and end  */
	public static String toBed3(final Locatable loc) {
		return String.join(
				"\t",
				loc.getContig(),
				String.valueOf(loc.getStart()-1) ,
				String.valueOf(loc.getEnd())
				);
		}
}
