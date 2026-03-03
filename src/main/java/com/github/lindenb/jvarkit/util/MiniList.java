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
package com.github.lindenb.jvarkit.util;


import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
/**
 * Minimal read-only List size() and get(index)
 * @param <T>
 */
public interface MiniList<T>  extends Iterable<T> {
	public default boolean isEmpty() {
		return size()==0;
		}
	public int size();
	public T get(int index);
	@Override
	default Iterator<T> iterator() {
		return new Iter<T>(this);
		}
	
	public default boolean contains(T o) {
		return indexOf(o)!=-1;
		}
	
	public default int indexOf(T o) {
		for(int i=0;i<size();i++) {
			if(Objects.equals(get(i), o)) return i;
			}
		return -1;
		}
	
	public default int lastIndexOf(T o) {
		for(int i=size()-1;i>=0;i--) {
			if(Objects.equals(get(i), o)) return i;
			}
		return -1;
		}
	 
    public default Stream<T> stream() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(),0), false);
    	}
	
	public default <F> MiniList<F> map(final Function<T,F> mapper) {
		final MiniList<T> delegate=this;
		return new MiniList<F>() {
			@Override
			public F get(int index) {
				return mapper.apply(delegate.get(index));
				}
			@Override
			public int size() {
				return delegate.size();
				}
			};
		}
	
	public default  MiniList<T> reverse() {
		final MiniList<T> delegate=this;
		return new MiniList<T>() {
			@Override
			public T get(int index) {
				return delegate.get((delegate.size()-1)-index);
				}
			@Override
			public int size() {
				return delegate.size();
				}
			};
		}
	
	
	static class Iter<E> implements Iterator<E>{
		private final MiniList<E> owner;
		private int idx=-1;
		Iter(MiniList<E> owner) {
			this.owner = owner;
			}
		@Override
		public boolean hasNext() {
			return idx+1 < owner.size();
			}
		@Override
		public E next() {
			if(this.idx+1>= owner.size()) throw new NoSuchElementException();
			this.idx++;
			return owner.get(idx);
			}
		}
	
	public static <T> MiniList<T> of(final List<T> delegate) {
		return new MiniList<T>() {
			@Override
			public T get(int index) {
				return delegate.get(index);
				}
			@Override
			public int size() {
				return delegate.size();
				}
			};
		}
	}
