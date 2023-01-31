/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.stream;

import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

public class CollectorsUtils {

	private static class MonoSet<T> extends AbstractSet<T> {
		private T value = null;
		@Override
		public boolean add(final T e) {
			if(e==null) throw new IllegalArgumentException("null value is not supported");
			if(this.value!=null ) {
				if(this.value==e || this.value.equals(e)) return false;
				throw new IllegalArgumentException("expected one unique value but got two:"+this.value+" "+e);
				}
			this.value=e;
			return true;
			}
		Optional<T> getOptional() {
			return Optional.ofNullable(this.value);
			}
		T getRequired() {
			if(this.value==null) throw new IllegalArgumentException("expected one  value but got none");
			return this.value;
			}
		@Override
		public Iterator<T> iterator() {
			return value==null?
					Collections.emptyIterator():
					Collections.singleton(value).iterator();
		}

		@Override
		public int size() {
			return (this.value==null?0:1);
		}
	}
	
/** expected one and only one value from this stream. Objects are compared using equals */
public static <T> Collector<T, ?, T> one() {
	return new DefaultCollector<>(
			()->new MonoSet<T>(),
			(SET,V)->SET.add(V),
			(A, B) -> { A.addAll(B); return A; },
			(SET)->SET.getRequired(),
			Collections.emptySet()
			);
	}
/** expected one or no value from this stream. Objects are compared using equals */
public static <T> Collector<T, ?, Optional<T>> optional() {
	return new DefaultCollector<>(
			()->new MonoSet<T>(),
			(SET,V)->SET.add(V),
			(A, B) -> { A.addAll(B); return A; },
			(SET)->SET.getOptional(),
			Collections.emptySet()
			);
	}
}
