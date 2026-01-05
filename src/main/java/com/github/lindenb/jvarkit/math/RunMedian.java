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
package com.github.lindenb.jvarkit.math;

import java.util.Arrays;

/** Compute running medians of odd span */
public class RunMedian {
private final int windowSize;

public RunMedian(final int windowSize) {
	if(windowSize<1 || windowSize%2!=1) throw new IllegalArgumentException("window size must be positive and odd");
	this.windowSize = windowSize;
	}

public int getWindowSize() {
	return windowSize;
	}

/** https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/runmed */
public static int getTurlachSize(final int n) {
	 return 1 + (int)(2 * Math.min((n-1)% 2, Math.ceil(0.1*n)));
}

public double[] applyToDoubleArray(final int[] input) {
	return apply(Arrays.stream(input).mapToDouble(I->(double)I).toArray());
	}

public int[] applyToIntArray(final int[] input) {
	return Arrays.stream(applyToDoubleArray(input)).mapToInt(I->(int)I).toArray();
	}


/** supid implementation.. */
public double[] apply(final double input[]) {
	final int w = getWindowSize();
	if(w==1 || input.length==0) {
		return Arrays.copyOf(input, input.length);
		}
	else
		{
		final Median median = new Median(w+1);
		final int halfw = w/2;
		final double[] dest = new double[input.length];
		for(int i=0;i< input.length;i++) {
			median.reset();
			final int beg_idx = Math.max(0,i-halfw);
			final int endg_idx = Math.min(i+halfw,input.length);
			for(int x=beg_idx;x< endg_idx;x++) {
				median.accept(input[i]);
				}
			dest[i] = median.getAsDouble();
			}
		
		return dest;
		}
	}

@Override
public String toString() {
	return "runmedian("+getWindowSize()+")";
	}
}
