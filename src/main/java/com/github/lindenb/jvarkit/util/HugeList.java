package com.github.lindenb.jvarkit.util;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public interface HugeList<T> extends Iterable<T>{
public long size();
public T get(long index);

public default HugeList<T> subList(final long fromIndex,final long toIndex) {
	final HugeList<T> me = this;
	return new HugeList<T>() {
		@Override public long size() { return toIndex - fromIndex;}
		@Override public T get(long index) { return me.get(fromIndex+index);}
		};
	}

public default Stream<T> stream() {
	return StreamSupport.stream(this.spliterator(), false);
	}

@Override
public default Iterator<T> iterator(){
	return new Iterator<T>()
		{
		long i = 0L;
		@Override
		public boolean hasNext() {
			return this.i< HugeList.this.size();
			}
		@Override
		public T next() {
			final T o = HugeList.this.get(this.i); 
			this.i++;
			return o;
			}
		};
	}

public static <T> HugeList<T> wrap(final List<T> list) {
	return new HugeList<T>() {
		@Override
		public T get(long index) {
			if(index>(long)Integer.MAX_VALUE)
				{
				throw new IllegalStateException();
				}
			return list.get((int)index);
			}
		@Override
		public long size() {
			return list.size();
			}
		};
	}

}
