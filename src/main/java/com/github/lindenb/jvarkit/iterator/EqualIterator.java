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
package com.github.lindenb.jvarkit.iterator;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.util.PeekableIterator;

public class EqualIterator<T> extends AbstractCloseableIterator<List<T>>
	{
	private final Comparator<T> comparator;
	private final PeekableIterator<T> delegate;
	public EqualIterator(final Iterator<T> delegate,final Comparator<T> comparator) {
		this.delegate = new PeekableIterator<>(delegate);
		this.comparator = comparator;
		}
	
	/** construct EqualIterator with default 'equals' comparator */
	public EqualIterator(final Iterator<T> delegate) {
		this(delegate,(A,B)->A.equals(B)?0:1);
		}
	
	protected boolean same(final T t1,final T t2) {
		return this.comparator.compare(t1, t2)==0;
	}
	
	@Override
	protected List<T> advance() {
		if(!delegate.hasNext()) return null;
		final List<T> L= new ArrayList<>();
		L.add(delegate.next());
		while(delegate.hasNext()) {
			final T rec2 = delegate.peek();
			if(!same(L.get(0),rec2)) {
				break;
			}
			L.add(delegate.next());
		}
		return L;
		}
	@Override
	public void close() {
		this.delegate.close();
		}
	}
