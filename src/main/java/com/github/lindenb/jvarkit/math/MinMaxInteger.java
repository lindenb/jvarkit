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
import java.util.OptionalInt;
import java.util.function.IntConsumer;
import java.util.function.IntUnaryOperator;
import java.util.stream.IntStream;

public class MinMaxInteger implements IntConsumer, IntUnaryOperator {
	private long count=0L;
	private int minV;
	private int maxV;
	
	/** allow this class to be used in stream, returns the original value */
	@Override
	public final int applyAsInt(int v) {
		accept(v);
		return v;
		}
	
	/** loop accept(v) over all array */
	public void accept(int...values) {
		for(int i=0;i< values.length;++i) {
			accept(values[i]);
			}
		}
	
	@Override
	public void accept(int value) {
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
	
	
	public MinMaxInteger() {
		}
	
	
	public MinMaxInteger(int m,int M) {
		if(m>M) throw new IllegalArgumentException("min "+m+" > max "+M);
		minV = m;
		maxV = M;
		count=2L;
		}

	public MinMaxInteger(int v) {
		minV = v;
		maxV = v;
		count=1L;
		}


	
	/** create a new instance with is the combination of 'this' and 'other' */
	public MinMaxInteger plus(final MinMaxInteger other) {
		
		if(this.count>0 && other.count>0) {
			final MinMaxInteger mM = new MinMaxInteger();
			mM.count = this.count+other.count;
			mM.minV = Math.min(this.minV,other.minV);
			mM.maxV = Math.max(this.maxV,other.maxV);
			return mM;
			}
		else if(this.count==0 && other.count==0) {
			return  new MinMaxInteger();
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
	
	public OptionalInt getMin() {
		return isEmpty()? OptionalInt.empty() : OptionalInt.of(minV);
		}
	
	public OptionalInt getMax() {
		return isEmpty()? OptionalInt.empty() : OptionalInt.of(maxV);
		}
	
	public int getMinAsInt() {
		if(isEmpty()) throw new IllegalStateException("no data was found");
		return minV;
		}
	
	public int getMaxAsInt() {
		if(isEmpty()) throw new IllegalStateException("no data was found");
		return maxV;
		}
	
	/** return a value between 0 and 1 */
	public double normalize(int value) {
		if(isEmpty()) throw new IllegalStateException("no data was found");
		return (value-this.minV)/(double)(this.maxV-this.minV); 
		}
	
	@Override
	public int hashCode() {
		if(this.count==0L)return 0;
		return Integer.hashCode(this.minV)*31 + Integer.hashCode(this.maxV);
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof MinMaxInteger)) return false;
		final MinMaxInteger mm = MinMaxInteger.class.cast(obj);
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
	protected MinMaxInteger clone() {
		final MinMaxInteger mM = new MinMaxInteger();
		mM.count=this.count;
		mM.minV=this.minV;
		mM.maxV=this.maxV;
		return mM;
		}
	
	public static MinMaxInteger of(final IntStream stream) {
		final MinMaxInteger mM = new MinMaxInteger();
		stream.forEach(mM);
		return mM;
		}
	
	public static MinMaxInteger of(int[] array,int offset, int length) {
		return of(Arrays.stream(array, offset,offset+length));
		}
	
	public static MinMaxInteger of(int[] array) {
		return of(array,0,array.length);
		}
	
	}
