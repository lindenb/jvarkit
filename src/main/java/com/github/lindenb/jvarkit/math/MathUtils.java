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

import java.math.BigInteger;

/**
 * 
 * utilities for Math
 *
 */
public class MathUtils {
	public static BigInteger factorial(long n) 
		{
		if(n<0L) throw new IllegalArgumentException("n<0");
	    BigInteger f = BigInteger.valueOf(1L);
	    for (long i = n; i > 0; i--) {
	          f = f.multiply(BigInteger.valueOf(i));
	    	}
	    return f;
		}

	/** https://en.wikipedia.org/wiki/Combination */
	public static BigInteger combination(long n,long k)
		{
		if(n-k<=0) throw new IllegalArgumentException("n-k<=0 with n="+n+" and k="+k);
		return factorial(n).divide(factorial(k).multiply(factorial(n-k)));
		}
	
	/** see R code for ppoints */
	public static  double[] ppoints(int n) {
		double a = n <= 10? 3.0/8 : 1.0/2;
		
		final double[] array=new double[n];
		for(int i=1;i<=n;i++) {
			array[i-1] = (i - a)/(n + 1.0 - 2*a);
			}				        
	    return array;
		}	
}
