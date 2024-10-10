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

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

/**
 * C++ like <algorithm>.
 * <code>List&lt;T&gt</code> using an element of type <code>List&lt;V&gt</code>
 *
 */
public class Algorithm<T,V> {
	
	private final Comparator<V> comparator;
	private final Function<T,V> extractor;
	
	@SuppressWarnings("unchecked")
	public Algorithm(final Function<T,V> extractor) {
		this(extractor,(A,B)->
			((Comparable<V>)A).compareTo(B)
			);
		}
	
	public Algorithm(final Function<T,V> extractor,final Comparator<V> comparator) {
		this.comparator = comparator;
		this.extractor = extractor;
		}

	private boolean lower_than(V a,V b) {
		return this.comparator.compare(a, b) < 0;
		}
	
	public int unsortedIndex(final List<T> dataVector) {
		if(dataVector.size()>1) {
			V prev = this.extractor.apply(dataVector.get(0));
			int i=1;
			while(i < dataVector.size()) {
				final V curr = this.extractor.apply(dataVector.get(i));
				if(this.lower_than(curr, prev)) {
					return i-1;
					}
				prev = curr;
				i++;
				}
			}
		return -1;
		}
	public boolean isSorted(final List<T> dataVector) {
		return unsortedIndex(dataVector)==-1;
		}
	
	public List<T>  assertSorted(final List<T> dataVector) {
		int i=unsortedIndex(dataVector);
		if(i!=-1){
			throw new IllegalArgumentException("input is not sorted at index="+i+"/"+dataVector.size());
			}
		return dataVector;
		}
	
	public void sort(final List<T> dataVector) {
		Collections.sort(dataVector, (A,B)->{
			return comparator.compare( extractor.apply(A), extractor.apply(B));
			});
		}
	public void sortIfNeeded(final List<T> dataVector) {
		if(!isSorted(dataVector)) {
			sort(dataVector);
			}
		}
	
	/** C+ lower_bound */
	public int lower_bound(
			final List<T> dataVector,
	        final V select
	        )
		{
		return lower_bound(dataVector,0,dataVector.size(), select);
		}

	/** C+ lower_bound */
	public int lower_bound(
				final List<T> dataVector,
				int first, 
				final int last,
				final V select
	            )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;

	            final T midObject  = dataVector.get(middle);
	            final V mid = this.extractor.apply(midObject);
	            if (this.lower_than(mid,select))
	                    {
	                    first = middle + 1;
	                    len = len - half - 1;
	                    }
	            else
	                    {
	                    len = half;
	                    }
	            }
	    return first;
	    }

	
	/** C+ upper_bound */
	public  int upper_bound(
			final List<T> dataVector,
	        final V select
	        )
		{
		return upper_bound(dataVector,0,dataVector.size(),select);
		}

	
	/** C+ upper_bound */
	public int upper_bound(
			final List<T> dataVector,
			int first,
			final int last,
			final V select
	        )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	    		final int half = len / 2;
	    		final int middle = first + half;

	            final T midObject = dataVector.get(middle);
	            final V mid = this.extractor.apply(midObject);
	            if (this.lower_than(select,mid))
	                    {
	                    len = half;
	                    }
	            else
	                    {
	                    first = middle + 1;
	                    len = len -  half - 1;
	                    }
	            }
	    return first;
	    }
	
	
	public int[] equal_range(
			final List<T> dataVector,
			final V select
			)
		{
		return equal_range(dataVector,0,dataVector.size(), select);
		}
	
	public int[] equal_range(
			final List<T> dataVector,
			int first,
			final int last,
			final V select
			) {
		int __len = last - first;

		while (__len > 0) {
			final int __half = __len / 2;
			int __middle = first + __half;
			final T mid = dataVector.get(__middle);
			final V midObj = this.extractor.apply(mid);
			if (this.lower_than(midObj, select)) {
				first = __middle;
				++first;
				__len = __len - __half - 1;
			} else if (this.lower_than(select,midObj)) {
				__len = __half;
			} else {
				int __left = lower_bound(dataVector, first, __middle, select);
				first += __len;
				++__middle;
				int __right = upper_bound(dataVector, __middle, first, select);
				return new int[] { __left, __right };
			}
		}
		return new int[] { first, first };
	}
	
	
	public List<T> equalList(
			final List<T> dataVector,
			final V lowerV,
			final V upperV
			)
		{
		return equalList(dataVector,0,dataVector.size(),lowerV,upperV);
		}

	
	public List<T> equalList(
			final List<T> dataVector,
			int first,
			final int last,
			final V lowerV,
			final V upperV
			)
		{
		if(lower_than(upperV, lowerV)) throw new IllegalArgumentException("upper < lowerV");
		final int i0 = this.lower_bound(dataVector, first, last, lowerV);
		if(i0<first) throw new IllegalStateException("index(lower)="+i0);
		if(i0>last) throw new IllegalStateException("index(lower)> vector.length="+dataVector.size());
		final int i1 = this.upper_bound(dataVector, i0, last, upperV);
		if(i1<first) throw new IllegalStateException("index(upper)="+i1);
		if(i1>last) throw new IllegalStateException("index(upper)> vector.length="+dataVector.size());
		if(i0>i1) throw new IllegalStateException("index(lower)="+i0+" > index(upper)="+i1);
		if(i0==i1) return Collections.emptyList();
		return dataVector.subList(i0,i1);
		}
	
	public List<T> equalList(
			final List<T> dataVector,
			final V select
			)
		{
		return equalList(dataVector,0,dataVector.size(),select);
		}
	
	public List<T> equalList(
			final List<T> dataVector,
			int first,
			final int last,
			final V select
			)
		{
		final int indexes[]= equal_range(dataVector, first, last, select);
		if(indexes[0]<0) throw new IllegalStateException("index[0]="+indexes[0]);
		if(indexes[1]<0) throw new IllegalStateException("index[1]="+indexes[1]);
		if(indexes[0]>indexes[1]) throw new IllegalStateException("index[1]="+indexes[1]+"<index[0]="+indexes[0]);
		if(indexes[0]==indexes[1]) return Collections.emptyList();
		return dataVector.subList(indexes[0], indexes[1]);
		}
	
	}
