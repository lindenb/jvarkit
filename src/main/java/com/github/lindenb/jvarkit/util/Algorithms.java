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
package com.github.lindenb.jvarkit.util;

import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import htsjdk.samtools.util.IterableAdapter;

/**
 * C++ like <algorithm>
 *
 */
public class Algorithms {
	
	private static <T >boolean lower_than(final T a, final T b,final Comparator<T> cmp)
		{
		return cmp.compare(a, b)<0;
		}
	
	/** C+ lower_bound */
	public static <T extends Comparable<T>> int lower_bound(
			final List<T> dataVector,
	        final T select
	        )
		{
		return lower_bound(dataVector,0,dataVector.size(),
				select,
				(A,B)->A.compareTo(B)
				);
		}

	/** C+ lower_bound */
	public  static <T> int lower_bound(
				final List<T> dataVector,
				final int first, 
				final int last,
				final T select,
				final Comparator<T> comparator
	            )
	    {
	    return lower_bound(dataVector,first,last,select,comparator,A->A);
	    }

	/** C+ lower_bound */
	public  static <T,U> int lower_bound(
				final List<T> dataVector,
				int first, 
				final int last,
				final U select,
				final Comparator<U> comparator,
				final Function<T,U> extractor
	            )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;

	            final T midObject  = dataVector.get(middle);
	            final U mid = extractor.apply(midObject);
	            if (lower_than(mid,select,comparator))
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
	public  static <T extends Comparable<T>> int upper_bound(
			final List<T> dataVector,
	        final T select
	        )
		{
		return upper_bound(dataVector,0,dataVector.size(), select,(A,B)->A.compareTo(B));
		}

	
	/** C+ upper_bound */
	public static <T> int upper_bound(
			final List<T> dataVector,
			int first,
			final int last,
			final T select,
			final Comparator<T> comparator
	        )
	    {
		return upper_bound(dataVector,first,last,select,comparator,A->A);
	    }

	/** C+ upper_bound */
	public static <T,U> int upper_bound(
			final List<T> dataVector,
			int first,
			final int last,
			final U select,
			final Comparator<U> comparator,
			final Function<T,U> extractor
	        )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	    		final int half = len / 2;
	    		final int middle = first + half;

	            final T midObject = dataVector.get(middle);
	            final U mid = extractor.apply(midObject);
	            if (lower_than(select,mid,comparator))
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
	
	
	public static <T extends Comparable<T>> int[] equal_range(
			final List<T> dataVector,
			final T select
			)
		{
		return equal_range(dataVector,0,dataVector.size(), select,(A,B)->A.compareTo(B),A->A);
		}
	
	public static <T,U> int[] equal_range(
			final List<T> dataVector,
			int first,
			final int last,
			final U select,
			final Comparator<U> comparator,
			final Function<T,U> converter
			) {
		int __len = last - first;

		while (__len > 0) {
			final int __half = __len / 2;
			int __middle = first + __half;
			final T mid = dataVector.get(__middle);
			final U midObj = converter.apply(mid);
			if (lower_than(midObj, select,comparator)) {
				first = __middle;
				++first;
				__len = __len - __half - 1;
			} else if (lower_than(select,midObj,comparator)) {
				__len = __half;
			} else {
				int __left = lower_bound(dataVector, first, __middle, select, comparator,converter);
				first += __len;
				++__middle;
				int __right = upper_bound(dataVector, __middle, first, select, comparator,converter);
				return new int[] { __left, __right };
			}
		}
		return new int[] { first, first };
	}
	
	public static <T,U> Stream<T> equal_range_stream(
			final List<T> dataVector,
			int first,
			final int last,
			final U select,
			final Comparator<U> comparator,
			final Function<T, U> converter
			)
		{
		final int indexes[]= equal_range(dataVector, first, last, select, comparator,converter);
		if(indexes[0]<0) throw new IllegalStateException("index[0]="+indexes[0]);
		if(indexes[1]<0) throw new IllegalStateException("index[1]="+indexes[1]);
		if(indexes[0]>indexes[1]) throw new IllegalStateException("index[1]="+indexes[1]+"<index[0]="+indexes[0]);
		if(indexes[0]==indexes[1]) return Stream.empty();
		return dataVector.subList(indexes[0], indexes[1]).stream();
		}
	
	public static <T,U> Iterator<T> equal_range_iterator(
			final List<T> dataVector,
			int first,
			final int last,
			final U select,
			final Comparator<U> comparator,
			final Function<T,U> converter
			)
		{
		return equal_range_stream(dataVector, first, last, select, comparator,converter).iterator();
		}
	
	
	// integer data ========================================================
	
	public static boolean isSorted(final int dataVector[]) {
		for(int i=0;i+1< dataVector.length;i++) {
			if(dataVector[i]>dataVector[i+1]) return false;
			}
		return true;
		}
	
	public static int upper_bound(final int dataVector[],final int select) {
		return upper_bound(dataVector, 0,dataVector.length,select);
		}
	
	public static int upper_bound(final int dataVector[],int first,int last,int select) {
		int len = last - first;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;

	            if (!(select < dataVector[middle]))
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
	
	public static int lower_bound(final int dataVector[],final int select) {
		return lower_bound(dataVector, 0,dataVector.length,select);
		}

	public static int lower_bound(final int dataVector[],int first,int last,final int select)
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;
	
	            if (dataVector[middle] < select)
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

	
}
