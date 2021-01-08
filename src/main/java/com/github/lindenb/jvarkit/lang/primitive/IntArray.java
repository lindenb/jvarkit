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
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.PrimitiveIterator;
import java.util.stream.IntStream;

/**
 * wrapper for an array of int.
 */
public class IntArray implements Cloneable,Serializable,Iterable<Integer> {
	private static final long serialVersionUID = 1L;
	private int mSize = 0;
	private int[] array;
	
	/** default constructor */
	public IntArray() {
		this(100);
	}
	/** copy constructor */
	public IntArray(final IntArray copy) {
		this(copy.array,0,copy.mSize);
	}

	/** constructor with defined buffer capacity */
	public IntArray(final int capacity) {
		if(this.mSize<0) throw new IllegalArgumentException("capacity <0 : "+ capacity);
		this.mSize = 0;
		this.array = new int[capacity];
	}
	
	/** create a copy of values from off with size=len */
	public IntArray(int values[],int off, int len) {
		this(len);
		System.arraycopy(values, off, this.array, 0, len);
		this.mSize = len;
		}
	/** create a copy of values */
	public IntArray(int values[]) {
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
	public IntArray slice(int off, int len) {
		if(off+len>size()) throw new IndexOutOfBoundsException("0<="+(off+len)+"<="+size());
		return new IntArray(this.array,off,len);
	}
	
	
	private int check(final int idx) {
		if(idx < 0 || idx >= this.mSize) {
			throw new IndexOutOfBoundsException("0<="+idx+"<idx="+this.mSize);
		}
		return idx;
	}
	
	/** set size to zero */
	public IntArray clear() {
		this.mSize = 0;
		return this;
	}

	/** get size */
	public int size() {
		return this.mSize;
	}
	
	/** return true if size==0 */
	public boolean isEmpty() {
		return size() == 0;
	}
	
	/** return index-th value */
	public int get(int index) {
		return this.array[check(index)];
	}

	/** set index-th value */
	public IntArray set(int index,int value) {
		this.array[check(index)]= value;
		return this;
	}
	
	private void ensure(int n) {
		final int avail = this.array.length - this.mSize;
		if(avail<= n) {
			final long L0x = (long)this.mSize+(long)n;
			final long L1x = Math.min((long)Integer.MAX_VALUE,(long)(this.array.length)+(long)(this.array.length)/2L);
			final long L2x = Math.max(L0x,L1x);
			if(L2x>=(long)Integer.MAX_VALUE) throw new IllegalStateException("Array too large "+ L2x);
			
			this.array  = Arrays.copyOf(this.array, (int)L2x);
			}
		}
	
	/** push bash the value */
	public int add(int value) {
		ensure(1);
		this.array[this.mSize] = value;
		this.mSize++;
		return value;
	}
	
	public IntArray addAll(final IntArray o) {
		return addAll(o.array,0,o.mSize);
		}
	
	public IntArray addAll(int[] values) {
		return addAll(values,0,values.length);
		}
	
	public IntArray addAll(final int[] values,int off,int len) {
		ensure(len);
		System.arraycopy(values, off, this.array, this.mSize, len);
		this.mSize+=len;
		return this;
		}
	
	public IntArray addAll(final Collection<Integer> col) {
		ensure(col.size());
		final Iterator<Integer> iter = col.iterator();
		while(iter.hasNext()) {
			this.array[this.mSize] = iter.next();
			this.mSize++;
			}
		return this;
	}
	// 0123456789
	public int remove(int idx) {
		int old = this.array[check(idx)];
		if(idx+1< this.mSize) {
			System.arraycopy(
				this.array, idx+1,
				this.array, idx,
				this.mSize - (idx+1));
			}
		this.mSize--;
		return old;
	}

	public int indexOf(int index,int value) {
		while(index < size()) {
			if(get(index)==value) return index;
			index++;
			}
		return -1;
		}
	public int indexOf(int value) {
		return isEmpty()?-1:indexOf(0,value);
		}
	
	public boolean contains(int value) {
		return indexOf(value)==-1;
	}
	
	public int insert(int index,int value) {
		ensure(1);
		if(index<= this.mSize) {
			System.arraycopy(
				this.array, index,
				this.array, index+1,
				this.mSize - index);
			}
		this.array[index]=value;
		this.mSize++;
		return value;
	}
	
 	
	public IntStream stream() {
		return Arrays.stream(this.array, 0, this.mSize);
	}
	
	@Override
	public PrimitiveIterator.OfInt iterator() {
		return stream().iterator();
	}
	
	/** sort this data */
	public IntArray sort() {
		Arrays.sort(this.array,0,this.mSize);
		return this;
	}
	
	/** convert to array. The array is a *copy* of the original data */
	public int[] toArray() {
		return Arrays.copyOf(this.array, this.mSize);
	}
	
	/** clone this object */
	public IntArray clone() {
		return new IntArray(this);
	}
	
	@Override
	public int hashCode() {
		int result = 0;
		for(int i=0;i< this.mSize;i++) {
	        result = 31 * result + this.array[i];
			}
		return result;
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
	
	public static IntArray read(DataInputStream in) throws IOException {
		final int n=in.readInt();
		IntArray vec = new IntArray(n);
		for(int i=0;i< n;i++) {
			vec.add(in.readInt());
			}
		return vec;
		}
	public static void write(DataOutputStream out,IntArray vec) throws IOException {
		out.writeInt(vec.size());
		for(int i=0;i< vec.size();i++) {
			out.writeInt(vec.get(i));
			}
		}
	}
