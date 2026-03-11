/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.util.CloseableIterator;

/**
 * Implementation of a chained iterator
 * @param <T>
 */
public class ChainedIterator<T> extends AbstractCloseableIterator<T> {
	private final LinkedList<Iterator<T>> delegates ;
	public ChainedIterator(final List<Iterator<T>> list) {
		this.delegates = new LinkedList<>(list);
		}
	
	@Override
	protected T advance() {
		while(!delegates.isEmpty()) {
			final Iterator<T> front = delegates.getFirst();
			if(!front.hasNext()) {
				close(front);
				delegates.pop();
				continue;
				}
			return front.next();
			}
		return null;
		}
	@Override
	public void close() {
		while(!delegates.isEmpty()) {
			close(delegates.pop());
			}
		}
	
	private void close(Iterator<T> t) {
		if(t!=null && (t instanceof CloseableIterator)) {
			CloseableIterator.class.cast(t).close();
			}
		}
}
