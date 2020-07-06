/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.TreeMap;

/*
 * see https://stackoverflow.com/questions/3052924
 * @author lindenb
 *
 */
public class DiscreteMedian<T extends Number>
	{
	private final TreeMap<T,Long> counter = new TreeMap<>();
	private long size = 0L;
	public DiscreteMedian() {
		}
	
	public long size() {
		return this.size;
		}
	
	public boolean isEmpty() {
		return counter.isEmpty();
		}
	
	public DiscreteMedian<T> clear() {
		this.size = 0L;
		this.counter.clear();
		return this;
	}
	
	public DiscreteMedian<T> add(final T value) {
		final long n = this.counter.getOrDefault(value, 0L);
		this.counter.put(value, n+1L);
		this.size++;
		return this;
		}
	
	public DiscreteMedian<T> add(final DiscreteMedian<T> other) {
		for(T t: other.counter.keySet()) {
			final long n1 = other.counter.get(t);
			final long n2 = this.counter.getOrDefault(t, 0L);
			this.counter.put(t, n1+n2);
			this.size+=n1;
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
					long prev_count = this.counter.getOrDefault(prev, 0L);
					this.counter.put(prev, prev_count+n);
					this.size+=n;
					}
				if(curr==null) break;
				prev=curr;
				n=0L;
				}
			n++;
			}
		
		return this;
		}
	
	public OptionalDouble getAverage() {
		if(isEmpty()) return OptionalDouble.empty();
		return OptionalDouble.of(
			counter.entrySet().
			stream().
			mapToDouble(KV->KV.getKey().doubleValue()* KV.getValue()).
			sum()/this.size
			);
		}
	
	public OptionalDouble getStandardDeviation() {
		final OptionalDouble avg = getAverage();
		if(!avg.isPresent()) return OptionalDouble.empty();
		final double v1 = avg.getAsDouble();
		double sum=0;
		for(Map.Entry<T,Long> key: this.counter.entrySet()) {
			final double v2 = Math.pow(key.getKey().doubleValue() - v1,2.0);
			for(long n1=0;n1< key.getValue();++n1) {
				sum+=v2;
				}
			}
		sum/=this.size;
		if(sum<=0) return OptionalDouble.empty();
		return OptionalDouble.of(Math.sqrt(sum));
	}
	
	/** get median for this dataset */
	public OptionalDouble getMedian() {
		if(isEmpty()) return OptionalDouble.empty();
		final long mid_x = this.size/2L;
		if (this.size%2==1) {
			long n=0;
			for(T key : this.counter.keySet())   {
				final long c = this.counter.get(key);
				if(mid_x< n+c) return OptionalDouble.of(key.doubleValue());
				n+=c;
				}
			}
		else {
			long n=0;
		 	final Iterator<T> iter = this.counter.keySet().iterator();
		 	while(iter.hasNext()) {
		 		T key =  iter.next();
		 		final long count = this.counter.get(key);
		 		if((mid_x-1)>= n+count)
		 			{
		 			n+=count;
		 			continue;
		 			}
		 		double v1= key.doubleValue();
		 		if((mid_x)< n+count) {
		 			return OptionalDouble.of(v1);
		 			}
		 		if(!iter.hasNext()) throw new IllegalStateException();
		 		key = iter.next();
				double v2= key.doubleValue();
				return OptionalDouble.of((v1+v2)/2.0);	
		 		}
		  }
		throw new IllegalStateException();
		}
	}
