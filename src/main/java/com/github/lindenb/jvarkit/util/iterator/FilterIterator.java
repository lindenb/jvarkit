/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.iterator;

import java.util.Iterator;
import java.util.function.Predicate;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class FilterIterator<E>
	extends AbstractIterator<E>
	implements CloseableIterator<E>
	{
	private final Iterator<E> delegate;
	private Predicate<E> predicate;
	public FilterIterator(
			final Iterator<E> delegate,
			final Predicate<E> predicate
			) {
		this.delegate = delegate;
		this.predicate = predicate;
	}
	
	@Override
	protected E advance() {
		while(this.delegate.hasNext())
			{
			final E item = this.delegate.next();
			if(this.predicate.test(item)) return item;
			}
		return null;
		}
	@Override
	public void close() {
		CloserUtil.close(this.delegate);
	}
	
}
