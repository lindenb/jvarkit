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
package com.github.lindenb.jvarkit.math.stats;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public abstract class Percentile {
	public enum Type {MIN,MAX,MEDIAN,AVERAGE,RANDOM,SUM}

protected Percentile() {
}

public double evaluate(double[] values, int start, int length)
	{
	return evaluate(Arrays.stream(values,start,start+length));
	}
public double evaluate(double[] values)
	{
	return evaluate(values,0,values.length);
	}

public double evaluate(int[] values, int start, int length)
	{
	return evaluate(Arrays.stream(values,start,start+length));
	}
public double evaluate(int[] values)
	{
	return evaluate(values,0,values.length);
	}

public abstract double evaluate(final DoubleStream stream);

public double evaluate(final IntStream stream) {
	return evaluate(stream.mapToDouble(I->I));
	}


public static Percentile median() { return _PMedian;}
public static Percentile average() { return _PAvg;}
public static Percentile max() { return _PMax;}
public static Percentile min() { return _PMin;}
public static Percentile random() { return _PRandom;}
public static Percentile sum() { return _PSum;}

public static Percentile of(final Type t) {
	switch(t) {
		case MIN: return min();
		case MAX: return max();
		case AVERAGE: return average();
		case MEDIAN: return median();
		case RANDOM: return random();
		case SUM: return sum();
		default: throw new IllegalStateException("bad Percentile type :"+t);
		}
	}

private static final  Percentile _PMin = new Percentile() {
	@Override
	public double evaluate(final DoubleStream stream) {
		return stream.min().getAsDouble();
		}
	
	@Override public String toString() {return "min";}
	};
	
private static final  Percentile _PMax = new Percentile() {
	@Override
	public double evaluate(final DoubleStream stream) {
		return stream.max().getAsDouble();
		}
	@Override public String toString() {return "max";}
	};

private static final  Percentile _PAvg = new Percentile() {
	@Override
	public double evaluate(final DoubleStream stream) {
		return stream.average().getAsDouble();
		}
	@Override public String toString() {return "average";}
	};

private static final  Percentile _PSum = new Percentile() {
	@Override
	public double evaluate(final DoubleStream stream) {
		return stream.sum();
		}
	@Override public String toString() {return "sum";}
	};
	
private static final  Percentile _PMedian = new Percentile() {
	@Override
	public double evaluate(final DoubleStream stream) {
		final double copy[]= stream.sorted().toArray();
		final int mid_x= copy.length/2;
		if(copy.length==1)
			{
			return copy[0];
			}
		else if(copy.length%2==0)
	        {
			return (copy[mid_x-1]+copy[mid_x])/2.0;
	        }
		else
	        {
	        return copy[mid_x];
	        }
		}
	@Override public String toString() {return "median";}
	};

private static final  Percentile _PRandom = new Percentile() {
	final Random rnd = new Random();
	
	private int index(int start, int length)
		{
		return start + this.rnd.nextInt(length);
		}
	public double evaluate(final DoubleStream stream)
		{
		return evaluate(stream.toArray());
		}
	@Override
	public double evaluate(double[] values, int start, int length) {
		return values[index(start,length)];
		}
	@Override
	public double evaluate(int[] values, int start, int length) {
		return values[index(start,length)];
		}
	@Override public String toString() {return "random";}
	};

}
