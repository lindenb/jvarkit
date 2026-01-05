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

public class ArrayResizer
	{
	private Percentile percentile = Percentile.MEDIAN;
	public enum Percentile {MEDIAN,AVERAGE}
	
	ArrayResizer percentile(final Percentile percentile) {
		this.percentile = percentile;
		return this;
		}
	
	public double[] resizeToDouble(double[] array,int newSize) {
		final double[] result = new double[newSize];
		if(newSize==0) {
			return result;
			}
		else if(newSize==array.length) {
			for(int i=0;i< result.length;i++) {
				result[i]=array[i];
				}
			return result;
			}
		else if(newSize < array.length) {
			final double f= (double)array.length/((double)newSize);//e.g array=1000,new=100, f=10
			for(int i=0;i< result.length;i++) {
				final int x0 = (int)((i  )*f);
				final int x1 = (int)((i+1)*f);
				result[i] = applyPercentile(array,x0,x1);
				}
			return result;
			}
		else {
			final double f= (double)array.length/((double)newSize);//e.g new=1000,array=100, f=0.1
			for(int i=0;i< result.length;i++) {
				final int x0 = (int)((i)*f);
				result[i] = array[x0];
				}
			return result;
			}
		}
	
	private double applyPercentile(double[] array,int i1,int i2) {
		switch(this.percentile) {
			case AVERAGE: 
				return percentileAverage(array, i1, i2);
			case MEDIAN:
			default:
				return percentileMedian(array,i1,i2);
			}
		
		}
	private double percentileAverage(double[] array,int i1,int i2) {
		double sum=0.0;
		int n=0;
		while(i1 < i2 && i1 < array.length) {
			if(Double.isNaN(array[i1])) continue;
			sum+=array[i1];
			i1++;
			n++;
			}
		if(n==0) return 0.0;
		return sum/n;
		}
	
	private double percentileMedian(double[] array,int i1,final int i2) {
		final double[] vec = new double[i2-i1];
		int n=0;
		while(i1 < i2 && i1 < array.length) {
			if(Double.isNaN(array[i1])) continue;
			vec[n] = array[i1];
			i1++;
			n++;
			}

		if(n==0) return 0.0;
		
		Arrays.sort(vec,0,n);
		final int mid = n/2;
		if(n%2==1 || mid==0) {
			return vec[mid];
			}
		else
			{
			return (vec[mid-1]+vec[mid])/2.0;
			}
		}
	}
