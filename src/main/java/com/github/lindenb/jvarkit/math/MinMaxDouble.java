/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.util.OptionalDouble;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.DoubleStream;

public class MinMaxDouble implements DoubleConsumer, DoubleUnaryOperator {
	private long count=0L;
	private double minV;
	private double maxV;
	
	/** allow this class to be used in stream, returns the original value */
	@Override
	public final double applyAsDouble(double v) {
		accept(v);
		return v;
		}
	
	@Override
	public void accept(double value) {
		if(Double.isNaN(value)) throw new IllegalArgumentException("cannot accept(NaN)");
		
		if(count==0L) {
			minV = value;
			maxV = value;
			}
		else
			{
			minV = Math.min(minV,value);
			maxV = Math.max(maxV,value);
			}
		count++;
		}
	
	
	public MinMaxDouble() {
		}
	
	public MinMaxDouble(double m,double M) {
		if(Double.isNaN(m)) throw new IllegalArgumentException("NaNa(m)");
		if(Double.isNaN(M)) throw new IllegalArgumentException("NaNa(M)");
		if(m>M) throw new IllegalArgumentException("min "+m+" > max "+M);
		minV = m;
		maxV = M;
		count=2L;
		}

	public MinMaxDouble(double v) {
		if(Double.isNaN(v)) throw new IllegalArgumentException("NaNa(v)");
		minV = v;
		maxV = v;
		count=1L;
		}

	
	/** create a new instance with is the combination of 'this' and 'other' */
	public MinMaxDouble plus(MinMaxDouble other) {
		
		if(this.count>0 && other.count>0) {
			final MinMaxDouble mM = new MinMaxDouble();
			mM.count = this.count+other.count;
			mM.minV = Math.min(this.minV,other.minV);
			mM.maxV = Math.max(this.maxV,other.maxV);
			return mM;
			}
		else if(this.count==0 && other.count==0) {
			return  new MinMaxDouble();
			}
		else if(this.count>0) {
			return this.clone();
			}
		else{
			return other.clone();
			}
		}
	
	/** return number of time this IntConsumer was visited */
	public long getCount() {
		return count;
		}
	
	/** return true if count==0 */
	public boolean isEmpty() {
		return getCount()==0;
		}
	
	public OptionalDouble getMin() {
		return isEmpty()? OptionalDouble.empty() : OptionalDouble.of(minV);
		}
	
	public OptionalDouble getMax() {
		return isEmpty()? OptionalDouble.empty() : OptionalDouble.of(maxV);
		}
	
	public double getMinAsDouble() {
		if(isEmpty()) throw new IllegalStateException("no data was found");
		return minV;
		}
	
	public double getMaxAsDouble() {
		if(isEmpty()) throw new IllegalStateException("no data was found");
		return maxV;
		}
	
	@Override
	public int hashCode() {
		if(this.count==0L)return 0;
		return Double.hashCode(this.minV)*31 + Double.hashCode(this.maxV);
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof MinMaxDouble)) return false;
		final MinMaxDouble mm = MinMaxDouble.class.cast(obj);
		return this.isEmpty()== mm.isEmpty() && 
				this.minV == mm.minV &&
				this.maxV == mm.maxV;
		}
	@Override
	public String toString() {
		if(isEmpty()) return "MinMax[undefined/undefined]";
		return "["+this.minV+"/"+this.maxV+"]";
		}
	
	@Override
	protected MinMaxDouble clone() {
		final MinMaxDouble mM = new MinMaxDouble();
		mM.count=this.count;
		mM.minV=this.minV;
		mM.maxV=this.maxV;
		return mM;
		}
	
	public static MinMaxDouble of(final DoubleStream stream) {
		final MinMaxDouble mM = new MinMaxDouble();
		stream.forEach(mM);
		return mM;
		}
	
	public static MinMaxDouble of(double[] array,int offset, int length) {
		return of(Arrays.stream(array, offset,offset+length));
		}
	
	public static MinMaxDouble of(double[] array) {
		return of(array,0,array.length);
		}
	
	}
