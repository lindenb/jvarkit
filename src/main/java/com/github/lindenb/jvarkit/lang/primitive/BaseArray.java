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


import java.io.Serializable;
import java.util.AbstractList;
import java.util.List;
import java.util.stream.Stream;

/**
 * wrapper for an array of primitives
 */
public abstract class BaseArray<T> implements Cloneable,Serializable,Iterable<T> {
	private static final long serialVersionUID = 1L;
	protected int mSize = 0;

	/** default constructor */
	protected BaseArray() {
		}

	/** get size */
	public final int size() {
		return this.mSize;
		}
	/** return true if size==0 */
	public final boolean isEmpty() {
		return size() == 0;
		}
	/** check array index */
	protected int check(final int idx) {
		if(idx < 0 || idx >= this.mSize) {
			throw new IndexOutOfBoundsException("0<="+idx+"<idx="+this.mSize);
			}
		return idx;
		}
	
	protected int extendSize(final int arrayLen,final int n) {
		// L0x : what is asked
		final long L0x = (long)this.mSize+(long)n;
		// L1x : arraylength * 1.5
		final long L1x = Math.min((long)Integer.MAX_VALUE,(long)(arrayLen)+(long)(arrayLen)/2L);
		// min of L0x L1x
		final long L2x = Math.max(L0x,L1x);
		if(L2x>=(long)Integer.MAX_VALUE) throw new IllegalStateException("Array too large "+ L2x);
		return (int)L2x;
		}
	
	
	protected abstract T getElementAt(int idx);
	protected abstract T setElementAt(int idx,final T value);
	protected abstract void addElement(final T value);
	protected abstract T removeElementAt(int idx);
	
	/** return a backed List for this data */ 
	public List<T> asList() {
		return new AbstractList<T>()
			{
			@Override
			public T remove(final int index)
				{
				return BaseArray.this.removeElementAt(index);
				}
			@Override
			public boolean add(final T e)
				{
				BaseArray.this.addElement(e);
				return true;
				}
			@Override
			public T set(final int index, final T element)
				{
				return BaseArray.this.setElementAt(index, element);
				}
			@Override
			public T get(final int index)
				{
				return BaseArray.this.getElementAt(index);
				}
			@Override
			public int size()
				{
				return BaseArray.this.size();
				}
			};
		}
		}
