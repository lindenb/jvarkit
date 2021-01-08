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

import java.util.Collections;
import java.util.Iterator;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;

/** simple abstract class over AbstractCloseableIterator and CloseableIterator */
public abstract class AbstractCloseableIterator<T> extends AbstractIterator<T> implements CloseableIterator<T> {

/** creates an empty closeable iterator */
public static <T> AbstractCloseableIterator<T> wrap(final Iterator<T> delegate, final Runnable onCloseOrNull) {
	return new AbstractCloseableIterator<T>() {
		@Override
		protected T advance() {
			return delegate.hasNext()?delegate.next():null;
			}
		@Override
		public void close() {
			if(onCloseOrNull!=null) onCloseOrNull.run();
			}
		};
}
	
/** creates an empty closeable iterator */
public static <T> AbstractCloseableIterator<T> empty() {
	return wrap(Collections.emptyIterator(),null); 
	}
}
