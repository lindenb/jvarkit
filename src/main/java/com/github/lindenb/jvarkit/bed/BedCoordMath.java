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
package com.github.lindenb.jvarkit.bed;

public class BedCoordMath {
	/** return distance between start and end */
	public static int getLength(final BedInterval bed) { 
		return getLength(bed.getBedStart(),bed.getBedEnd());
	}
	/** return distance between start and end */
	public static int getLength(int start0, int end0) { 
		return end0 - start0;
	}
	
	/** tell whether two BED intervals overlap  */
    public static boolean overlaps(
    		final BedInterval r1,
    		final BedInterval r2
    		) {
    	return overlaps(
    			r1.getBedStart(),r1.getBedEnd(),
    			r2.getBedStart(),r2.getBedEnd()
    			) && r1.getContig().equals(r2.getContig());
    	}
	
	/** tell whether two BED intervals overlap  */
    public static boolean overlaps(
    		final int start, final int end,
    		final int start2, final int end2
    		) {
    	return !((end2 <= start) || (end<= start2));
    	}
    
    
    /**
     * Determines the amount of overlap between two coordinate ranges. Assumes that the two ranges
     * actually do overlap and therefore may produce strange results when they do not!
     */
    public static int getOverlap(final int start, final int end, final int start2, final int end2) {
        return getLength(Math.max(start, start2), Math.min(end, end2));
    }

}
