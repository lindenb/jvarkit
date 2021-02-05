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
package com.github.lindenb.jvarkit.lang.primitive;


import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.PrimitiveIterator;
import java.util.stream.DoubleStream;

/**
 * wrapper for an array of double.
 * GENERATED : DO NOT EDIT.
 */
public class DoubleArray extends BaseArray<Double> {
	private static final long serialVersionUID = 1L;
	private double[] array;
	
	/** default constructor */
	public DoubleArray() {
		this(100);
	}
	/** copy constructor */
	public DoubleArray(final DoubleArray copy) {
		this(copy.array,0,copy.mSize);
	}

	/** constructor with prefilled 'N' values  */
	public DoubleArray(final int size,final double defaultValue) {
		this(size);
		super.mSize = size;
		Arrays.fill(this.array,defaultValue);
	}


	/** constructor with defined buffer capacity */
	public DoubleArray(final int capacity) {
		super();
		if(capacity < 0) throw new IllegalArgumentException("capacity <0 : "+ capacity);
		this.array = new double[capacity];
	}
	
	/** create a copy of values from off with size=len */
	public DoubleArray(final double values[],int off, int len) {
		this(len);
		System.arraycopy(values, off, this.array, 0, len);
		super.mSize = len;
		}
	/** create a copy of values */
	public DoubleArray(final double values[]) {
		this(values,0,values.length);
	}
	
	/** create a IntArray by wrapping the already existing array 'values'
	 * returned IntArray is now owner of the array
	 * */
	public static DoubleArray wrap(final double values[]) {
		final DoubleArray a = new DoubleArray(0);
		a.array = values;
		a.mSize = values.length;
		return a;
	}
	
	/** slice a copy of the array */
	public DoubleArray slice(final int off,final int len) {
		if(off+len>size()) throw new IndexOutOfBoundsException("0<="+(off+len)+"<="+size());
		return new DoubleArray(this.array,off,len);
	}
	
	/** slice a copy of the array.  It extracts through the end of the sequence  */
	public DoubleArray slice(final int off) {
		return new DoubleArray(this.array,off,size()-off);
	}
	
	/** set size to zero */
	public DoubleArray clear() {
		super.mSize = 0;
		return this;
	}
	
	/** return index-th value */
	public double get(int index) {
		return this.array[check(index)];
	}

	/** set index-th value , return previous value*/
	public double set(int index,final double value) {
		final double old = this.array[check(index)];
		this.array[index] = value;
		return old;
	}
	
	private void ensure(final int n) {
		final int avail = this.array.length - super.mSize;
		if(avail < n) {
			this.array  = Arrays.copyOf(this.array,extendSize(this.array.length,n));
			}
		}
	
	/** push back the value */
	public double add(final double value) {
		ensure(1);
		this.array[super.mSize] = value;
		super.mSize++;
		return value;
	}
	
	public DoubleArray addAll(final DoubleArray o) {
		return addAll(o.array,0,o.mSize);
		}
	
	public DoubleArray addAll(final double[] values) {
		return addAll(values,0,values.length);
		}
	
	public DoubleArray addAll(final double[] values,int off,int len) {
		ensure(len);
		System.arraycopy(values, off, this.array, super.mSize, len);
		super.mSize+=len;
		return this;
		}
	
	public DoubleArray addAll(final Collection<Double> col) {
		ensure(col.size());
		final Iterator<Double> iter = col.iterator();
		while(iter.hasNext()) {
			this.array[super.mSize] = iter.next();
			super.mSize++;
			}
		return this;
	}
	
	/** remove 1st value */
	public double popFront() {
		return remove(0);
		}
	
	/** remove last value */
	public double popBack() {
		return remove(size()-1);
		}
	
	
	/** remove idx-th value */
	public double remove(final int idx) {
		final double old = this.array[check(idx)];
		if(idx+1< super.mSize) {
			System.arraycopy(
				this.array, idx+1,
				this.array, idx,
				super.mSize - (idx+1));
			}
		super.mSize--;
		return old;
	}

	public int indexOf(int index,final double value) {
		while(index < size()) {
			if(get(index)==value) return index;
			index++;
			}
		return -1;
		}
		
	public int indexOf(final double value) {
		return isEmpty()?-1:indexOf(0,value);
		}
	
	public boolean contains(final double value) {
		return indexOf(value)!=-1;
	}
	
	public double insert(int index,double value) {
		ensure(1);
		if(index<= super.mSize) {
			System.arraycopy(
				this.array, index,
				this.array, index+1,
				super.mSize - index);
			}
		this.array[index]=value;
		super.mSize++;
		return value;
	}
	
		public DoubleStream stream() {
		return Arrays.stream(this.array, 0, super.mSize);
	}
	
	
	@Override
	public PrimitiveIterator.OfDouble iterator() {
		return stream().iterator();
	}
		
	/** sort this data */
	public DoubleArray sort() {
		Arrays.sort(this.array,0,super.mSize);
		return this;
	}
	
	/** convert to array. The array is a *copy* of the original data */
	public double[] toArray() {
		return Arrays.copyOf(this.array, super.mSize);
	}
	
	/** clone this object */
	public DoubleArray clone() {
		return new DoubleArray(this);
	}
	
	@Override
	public int hashCode() {
		int result = 0;
		for(int i=0;i< super.mSize;i++) {
	        result = 31 * result + Double.hashCode(this.array[i]);
			}
		return result;
		}
		
	@Override
	protected final Double getElementAt(final int idx)
		{
		return get(idx);
		}
	@Override
	protected final Double setElementAt(final int idx,final Double value)
		{
		return set(idx,value);
		}
	@Override
	protected final void addElement(final Double value)
		{
		this.add(value);
		}
	@Override
	protected final Double removeElementAt(final int idx)
		{
		return this.remove(idx);
		}
	
	
	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		for(int i=0;i< size();i++) {
			if(i>0) sb.append(',');
			sb.append(get(i));
			}
		return sb.toString();
		}
	
	public static DoubleArray read(final DataInputStream in) throws IOException {
		final int n=in.readInt();
		final DoubleArray vec = new DoubleArray(n);
		for(int i=0;i< n;i++) {
			vec.add(in.readDouble());
			}
		return vec;
		}
	
	public static void write(final DataOutputStream out,DoubleArray vec) throws IOException {
		out.writeInt(vec.size());
		for(int i=0;i< vec.size();i++) {
			out.writeDouble(vec.get(i));
			}
		}
	}
