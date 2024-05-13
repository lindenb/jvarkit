/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util;

import java.util.AbstractList;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Like a list with long index support
 *
 * @param <T>
 */
public interface LongList<T> extends Iterable<T> {
	public long size();
	public default boolean isEmpty() {
		return size()==0;
		}
	public T get(long idx);
	public default T remove(long idx) {
		throw new UnsupportedOperationException();
		}
	
	public default Iterator<T> iterator() {
		return new Iterator<T>() {
			long i=0L;
			@Override
			public boolean hasNext() {
				return i < size();
				}
			@Override
			public T next() {
				final T o = get(i);
				i++;
				return o;
				}
			};
		}
	public default Stream<T> stream() {
        final Spliterator<T> s = Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED);
        return StreamSupport.stream(s, false);
		}
	
	public static <T> LongList<T> of(final List<T> list) {
		return new LongList<T>() {
			@Override
			public T get(long idx) {
				if(idx<0 || idx>=(long)list.size()) throw new IndexOutOfBoundsException();
				return list.get((int)idx);
				}
			@Override
			public long size() {
				return  list.size();
				}
			};
		}
	
	public default List<T> asList() {
		if(size()>= Integer.MAX_VALUE) throw new IndexOutOfBoundsException("too many items to convert to List");
		return new AbstractList<T>() {
			@Override
			public int size() {
				return (int)LongList.this.size();
				}
			@Override
			public T get(int index) {
				return LongList.this.get((long)index);
				}
			};
		}
	}
