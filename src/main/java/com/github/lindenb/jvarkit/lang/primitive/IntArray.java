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
import java.util.stream.IntStream;

/**
 * wrapper for an array of int.
 * GENERATED : DO NOT EDIT.
 */
public class IntArray extends BaseArray<Integer> {
	private static final long serialVersionUID = 1L;
	private int[] array;
	
	/** default constructor */
	public IntArray() {
		this(100);
	}
	/** copy constructor */
	public IntArray(final IntArray copy) {
		this(copy.array,0,copy.mSize);
	}

	/** constructor with prefilled 'N' values  */
	public IntArray(final int size,final int defaultValue) {
		this(size);
		super.mSize = size;
		Arrays.fill(this.array,defaultValue);
	}


	/** constructor with defined buffer capacity */
	public IntArray(final int capacity) {
		super();
		if(capacity < 0) throw new IllegalArgumentException("capacity <0 : "+ capacity);
		this.array = new int[capacity];
	}
	
	/** create a copy of values from off with size=len */
	public IntArray(final int values[],int off, int len) {
		this(len);
		System.arraycopy(values, off, this.array, 0, len);
		super.mSize = len;
		}
	/** create a copy of values */
	public IntArray(final int values[]) {
		this(values,0,values.length);
	}
	
	/** create a IntArray by wrapping the already existing array 'values'
	 * returned IntArray is now owner of the array
	 * */
	public static IntArray wrap(final int values[]) {
		final IntArray a = new IntArray(0);
		a.array = values;
		a.mSize = values.length;
		return a;
	}
	
	/** slice a copy of the array */
	public IntArray slice(final int off,final int len) {
		if(off+len>size()) throw new IndexOutOfBoundsException("0<="+(off+len)+"<="+size());
		return new IntArray(this.array,off,len);
	}
	
	/** slice a copy of the array.  It extracts through the end of the sequence  */
	public IntArray slice(final int off) {
		return new IntArray(this.array,off,size()-off);
	}
	
	/** set size to zero */
	public IntArray clear() {
		super.mSize = 0;
		return this;
	}
	
	/** return index-th value */
	public int get(int index) {
		return this.array[check(index)];
	}

	/** set index-th value , return previous value*/
	public int set(int index,final int value) {
		final int old = this.array[check(index)];
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
	public int add(final int value) {
		ensure(1);
		this.array[super.mSize] = value;
		super.mSize++;
		return value;
	}
	
	public IntArray addAll(final IntArray o) {
		return addAll(o.array,0,o.mSize);
		}
	
	public IntArray addAll(final int[] values) {
		return addAll(values,0,values.length);
		}
	
	public IntArray addAll(final int[] values,int off,int len) {
		ensure(len);
		System.arraycopy(values, off, this.array, super.mSize, len);
		super.mSize+=len;
		return this;
		}
	
	public IntArray addAll(final Collection<Integer> col) {
		ensure(col.size());
		final Iterator<Integer> iter = col.iterator();
		while(iter.hasNext()) {
			this.array[super.mSize] = iter.next();
			super.mSize++;
			}
		return this;
	}
	
	/** remove 1st value */
	public int popFront() {
		return remove(0);
		}
	
	/** remove last value */
	public int popBack() {
		return remove(size()-1);
		}
	
	
	/** remove idx-th value */
	public int remove(final int idx) {
		final int old = this.array[check(idx)];
		if(idx+1< super.mSize) {
			System.arraycopy(
				this.array, idx+1,
				this.array, idx,
				super.mSize - (idx+1));
			}
		super.mSize--;
		return old;
	}

	public int indexOf(int index,final int value) {
		while(index < size()) {
			if(get(index)==value) return index;
			index++;
			}
		return -1;
		}
		
	public int indexOf(final int value) {
		return isEmpty()?-1:indexOf(0,value);
		}
	
	public boolean contains(final int value) {
		return indexOf(value)!=-1;
	}
	
	public int insert(int index,int value) {
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
	
		public IntStream stream() {
		return Arrays.stream(this.array, 0, super.mSize);
	}
	
	
	@Override
	public PrimitiveIterator.OfInt iterator() {
		return stream().iterator();
	}
		
	/** sort this data */
	public IntArray sort() {
		Arrays.sort(this.array,0,super.mSize);
		return this;
	}
	
	/** convert to array. The array is a *copy* of the original data */
	public int[] toArray() {
		return Arrays.copyOf(this.array, super.mSize);
	}
	
	/** clone this object */
	public IntArray clone() {
		return new IntArray(this);
	}
	
	@Override
	public int hashCode() {
		int result = 0;
		for(int i=0;i< super.mSize;i++) {
	        result = 31 * result + Integer.hashCode(this.array[i]);
			}
		return result;
		}
		
	@Override
	protected final Integer getElementAt(final int idx)
		{
		return get(idx);
		}
	@Override
	protected final Integer setElementAt(final int idx,final Integer value)
		{
		return set(idx,value);
		}
	@Override
	protected final void addElement(final Integer value)
		{
		this.add(value);
		}
	@Override
	protected final Integer removeElementAt(final int idx)
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
	
	public static IntArray read(final DataInputStream in) throws IOException {
		final int n=in.readInt();
		final IntArray vec = new IntArray(n);
		for(int i=0;i< n;i++) {
			vec.add(in.readInt());
			}
		return vec;
		}
	
	public static void write(final DataOutputStream out,IntArray vec) throws IOException {
		out.writeInt(vec.size());
		for(int i=0;i< vec.size();i++) {
			out.writeInt(vec.get(i));
			}
		}
	}
