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
package com.github.lindenb.jvarkit.math;

import java.util.OptionalDouble;
import java.util.function.DoubleConsumer;
import java.util.function.Supplier;

/**
 * Average: stupid wrapper around count/sum
 */
public class Average implements Supplier<OptionalDouble>,DoubleConsumer,Comparable<Average> {
	private long count = 0L;
	private double sum = 0.0;
	public Average() {
		}
	/** copy constructor */
	public Average(final Average cp) {
		this.count = cp.count;
		this.sum = cp.sum;
		}
	
	public void reset() {
		this.count = 0L;
		this.sum = 0.0;
		}
	
	public long getCount() {
		return this.count;
		}
	public double getTotal() {
		return this.sum;
		}
	
	@Override
	public int compareTo(final Average o) {
		final OptionalDouble a =  get();
		final OptionalDouble b= o.get();
		if(a.isPresent() && !b.isPresent()) return -1;
		if(!a.isPresent() && b.isPresent()) return  1;
		if(!a.isPresent() && !b.isPresent()) return  0;
		return Double.compare(a.getAsDouble(), b.getAsDouble());
		}
	
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Average)) return false;
		return this.compareTo(Average.class.cast(obj))==0;
		}
	
	public double getAsDouble() {
		if(this.count==0L) throw new IllegalStateException("count==0");
		return sum/count;
		}
	
	@Override
	public int hashCode() {
		return get().hashCode();
		}
	
	@Override
	public void accept(final double value) {
		sum+=value;
		count++;
		}
	@Override
	public OptionalDouble get() {
		return count==0L?
			OptionalDouble.empty():
			OptionalDouble.of(sum/count)
			;
		}
	
	@Override
	protected Average clone() {
		return new Average(this);
		}
	
	@Override
	public String toString() {
		OptionalDouble v = get();
		return v.isPresent()?String.valueOf(v.getAsDouble()):"N/A";
		}
	}
