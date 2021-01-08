/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.DoubleStream;
/*
 * see https://stackoverflow.com/questions/3052924
 * @author lindenb
 *
 */
public class DiscreteMedian<T extends Number>
	{
	public static enum Tendency  {
		median,
		average,
		min,
		max
		};
		
	private class Entry {
		final T key;
		long count = 0L;
		Entry(final T key) { this.key = key;}
		}
	private final List<Entry> counter = new ArrayList<>(1000);
	public DiscreteMedian() {
		}
	
	@Override
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append("size:").append(size()).append(" ");
		OptionalDouble opt = getMedian();
		sb.append("median:");
		if(opt.isPresent()) {
			sb.append(opt.getAsDouble());
			}
		else
			{
			sb.append("N/A");
			}
		opt = getAverage();
		sb.append(" average :");
		if(opt.isPresent()) {
			sb.append(opt.getAsDouble());
			}
		else
			{
			sb.append("N/A");
			}
		sb.append("\n");
		for(final Entry kv : this.counter) {
			sb.append("\t[").append(kv.key).append("]\t").append(kv.count).append("\n");
			}
		sb.append("\n");
		return sb.toString();
		}
	public long size() {
		return  this.counter.stream().mapToLong(KV->KV.count).sum();
		}
	
	public boolean isEmpty() {
		return counter.isEmpty();
		}
	
	public DiscreteMedian<T> clear() {
		this.counter.clear();
		return this;
	}
	
	@SuppressWarnings("unchecked")
	private int lower_bound(final T select) {
	    int len = this.counter.size();
	    int first = 0;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;

	            final T mid  = this.counter.get(middle).key;
	            if (Comparable.class.cast(mid).compareTo(select)<0)
	                    {
	                    first = middle + 1;
	                    len = len - half - 1;
	                    }
	            else
	                    {
	                    len = half;
	                    }
	            }
	    return first;
		}
	
	private void add(final T value,long incr) {
		final int idx = lower_bound(value);
		final Entry entry;
		if(idx>=this.counter.size() || !counter.get(idx).key.equals(value)) {
			entry = new Entry(value);
			this.counter.add(idx, entry);
			}
		else
			{
			entry = this.counter.get(idx);
			}
		entry.count+=incr;
		}
	
	public DiscreteMedian<T> add(final T value) {
		add(value,1L);
		return this;
		}
	
	public DiscreteMedian<T> add(final DiscreteMedian<T> other) {
		for(final Entry t: other.counter) {
			this.add(t.key,t.count);
			}
		return this;
		}
	
	public DiscreteMedian<T> addAll(final Collection<T> col) {
		Iterator<T> iter = col.iterator();
		T prev=null;
		long n=0L;
		for(;;) {
			T curr = iter.hasNext()?iter.next():null; 
			if(curr==null || (prev!=null && !prev.equals(curr))) {
				if(prev!=null) {
					this.add(prev,n);
					}
				if(curr==null) break;
				prev=curr;
				n=0L;
				}
			n++;
			}
		
		return this;
		}
	
	public OptionalDouble getTendency(final Tendency tendency) {
		switch(tendency) {
		case average: return getAverage();
		case max: return getMax();
		case median: return getMedian();
		case min: return getMin();
		default: throw new IllegalStateException();
		}
	}
	
	public DoubleStream getAsDoubleStream() {
		return  this.counter.stream().flatMapToDouble(
				KV->DoubleStream.generate(()->KV.key.doubleValue()).limit(KV.count)
				);
		}
	
	public OptionalDouble getMin() {
		return getAsDoubleStream().min();
		}
	
	public OptionalDouble getMax() {
		return getAsDoubleStream().max();
		}
	
	public OptionalDouble getAverage() {
		return getAsDoubleStream().average();
		}
	
	public OptionalDouble getStandardDeviation() {
		final OptionalDouble avg = getAverage();
		if(!avg.isPresent()) return OptionalDouble.empty();
		final double v1 = avg.getAsDouble();
		
		final OptionalDouble optdbl =  getAsDoubleStream().map(V-> Math.pow(V- v1,2.0)).average();
		if(!optdbl.isPresent() || optdbl.orElse(0.0)<=0.0) return OptionalDouble.empty();
		return OptionalDouble.of(Math.sqrt(optdbl.getAsDouble()));
		}
	
	/** get median for this dataset */
	public OptionalDouble getMedian() {
		if(isEmpty()) return OptionalDouble.empty();
		final long L= this.size();
		final long mid_x = L/2L;
		if (L%2==1) {
			// if size==1 , midx_x=0, we peek mid_x
			// if size==3 , midx_x=1, we peek mid_x
			return getAsDoubleStream().skip(mid_x).findFirst();
			}
		else {
			// if size==2 , midx_x=1, we peek 0 and 1
			// if size==4, mid_dx=2, we peek 1 and 2
			return getAsDoubleStream().skip(mid_x-1L).limit(2).average();
		  }
		}
	}
