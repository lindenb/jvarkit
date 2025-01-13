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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;

import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.function.DoubleConsumer;
import java.util.function.Supplier;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * Median calculation
 */
public class Median implements Supplier<OptionalDouble>,DoubleConsumer,Comparable<Median> {
	private static final int SIZE_SORT_ON_DISK= 2_000_000;
	private int count = 0;
	private double[] array;
	private boolean sorted=true;
	public Median() {
		this(10);
		}
	public Median(int buffer_size) {
		this.array = new double[buffer_size];
		this.count = 0;
		this.sorted = true;
		}
	/** copy constructor */
	public Median(final Median cp) {
		this.count = cp.count;
		this.array =  Arrays.copyOf(cp.array, cp.count);
		this.sorted = cp.sorted;
		}
	
	public void reset() {
		this.count = 0;
		Arrays.fill(this.array, 0.0);
		this.sorted = true;
		}
	
	public long getCount() {
		return this.count;
		}
	
	public OptionalDouble getAverage() {
		return Arrays.stream(this.array,0,this.count).average();
		}
	
	@Override
	public int compareTo(final Median o) {
		final OptionalDouble a =  get();
		final OptionalDouble b= o.get();
		if(a.isPresent() && !b.isPresent()) return -1;
		if(!a.isPresent() && b.isPresent()) return  1;
		if(!a.isPresent() && !b.isPresent()) return  0;
		return Double.compare(a.getAsDouble(), b.getAsDouble());
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Median)) return false;
		return this.compareTo(Median.class.cast(obj))==0;
		}
	
	public double getAsDouble() {
		return get().orElseThrow(()->new IllegalStateException("count==0"));
		}
	
	@Override
	public int hashCode() {
		return get().hashCode();
		}
	
	private void ensureSize(int new_size) {
		if(new_size >= array.length) {
			final long new_sizeL = Math.max((long)(array.length*1.5),(long)new_size);
			if(new_sizeL>=Integer.MAX_VALUE) throw new IllegalStateException("cannot resize array to "+new_sizeL);
			int new_size2  = (int)new_sizeL;
			this.array = Arrays.copyOf(this.array, new_size2);
			}
		}
	public void accept(final double[] values) {
		accept(values,0,values.length);
		}
	public void accept(double[] values,int start,int length) {
		if(length==0) return;
		ensureSize(this.count+length);
		for(int i=0;i< length;i++) {
			this.array[count++]=values[start+i];
			}
		sorted=false;
		}
	
	
	@Override
	public void accept(double value) {
		ensureSize(this.count+1);
		this.array[count]=value;
		sorted = false;
		count++;
		}
	
	private void sort() {
		if(!sorted) {
			if(this.count < SIZE_SORT_ON_DISK ) {
				Arrays.sort(this.array,0,this.count);
				}
			else
				{
				/** SortingLongCollection : cannot use this because the ordering of the double packed as long won't be good */
				class DoubleCodec extends AbstractDataCodec<Double> {
					@Override
					public Double decode(DataInputStream dis) throws IOException {
						try {
							return dis.readDouble();
							}
						catch(EOFException err) {
							return null;
							}
						catch(Throwable err) {
							throw new RuntimeIOException(err);
							}
						}
					@Override
					public void encode(DataOutputStream dos, Double object) throws IOException {
						dos.writeDouble(object);
						}
					@Override
					public DoubleCodec clone() {
						return new DoubleCodec();
						}
					}
				final SortingCollection<Double> sorter = SortingCollection.newInstance(
						Double.class,
						new DoubleCodec(),
						(A,B)->Double.compare(A,B),
						SIZE_SORT_ON_DISK,
						IOUtils.getDefaultTempDir()
						);
				sorter.setDestructiveIteration(true);
				for(int i=0;i< this.count;i++) {
					sorter.add(this.array[i]);
					}
				sorter.doneAdding();
				try(CloseableIterator<Double> iter=sorter.iterator()) {
					int i=0;
					while(iter.hasNext()) {
						this.array[i++] =  iter.next();
						}
					}
				sorter.cleanup();
				}
			sorted=true;
		}
	}
	
	@Override
	public OptionalDouble get() {
		if(count==0) return OptionalDouble.empty();
		sort();
		if(this.count==1)
			{
			return OptionalDouble.of(this.array[0]);
			}
		final int mid_x= this.count/2;
		if(this.count%2==0)
	        {
			return OptionalDouble.of((this.array[mid_x-1]+this.array[mid_x])/2.0);
	        }
		else
	        {
	        return OptionalDouble.of(this.array[mid_x]);
	        }
		}
	@Override
	public String toString() {
		final OptionalDouble v = get();
		return v.isPresent()?String.valueOf(v.getAsDouble()):"N/A";
		}
	}
