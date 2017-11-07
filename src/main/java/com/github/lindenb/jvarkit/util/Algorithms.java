/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import htsjdk.samtools.util.IterableAdapter;

/**
 * C++ like <algorithm>
 *
 */
public class Algorithms {
	/** C+ lower_bound */
	public static <T extends Comparable<T>> int lower_bound(
			List<T> dataVector,
	        T select
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
				int first, 
				final int last,
				final T select,
				final Comparator<T> comparator
	            )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	            final int half = len / 2;
	            final int middle = first + half;
	            final T x= dataVector.get(middle);
	            if (comparator.compare(x,select)<0)
	                    {
	                    first = middle + 1;
	                    len -= half + 1;
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
	public static <T extends Comparable<T>> int upper_bound(
			final List<T> dataVector,
			int first,
			final int last,
			final T select,
			final Comparator<T> comparator
	        )
	    {
	    int len = last - first;
	    while (len > 0)
	            {
	    		final int half = len / 2;
	    		final int middle = first + half;
	            final T x= dataVector.get(middle);
	            if (comparator.compare(select,x)<0)
	                    {
	                    len = half;
	                    }
	            else
	                    {
	                    first = middle + 1;
	                    len -= half + 1;
	                    }
	            }
	    return first;
	    }
	
	
	public static <T extends Comparable<T>> int[] equal_range(
			final List<T> dataVector,
			final T select
			)
		{
		return equal_range(dataVector,0,dataVector.size(), select,(A,B)->A.compareTo(B));
		}
	
	public static <T extends Comparable<T>> int[] equal_range(
			final List<T> dataVector,
			int first,
			final int last,
			final T select,
			final Comparator<T> comparator
			)
		{
		int __len = last - first;

	        while (__len > 0)
		  	{
		  	  final int __half = __len / 2;
		  	  int __middle = first + __half;
		  	  if (comparator.compare(dataVector.get(__middle), select)<0)
		  	    {
		  	      first = __middle;
		  	      ++first;
		  	      __len = __len - __half - 1;
		  	    }
		  	  else if (comparator.compare(select, dataVector.get(__middle))<0)
		  	  	{
		  	    __len = __half;
		  	  	}
		  	  else
		  	    {
		  	    int __left = lower_bound(dataVector,first, __middle, select, comparator);
		  	    first += __len;
		  	    ++__middle;
		  	    int __right = upper_bound(dataVector,__middle, first, select, comparator);
		  	    return new int[] {__left, __right};
		  	    }
		  	}
	    return new int[] {first, first};
		}
	
	public static <T extends Comparable<T>> Iterator<T> equal_range_iterator(
			final List<T> dataVector,
			int first,
			final int last,
			final T select,
			final Comparator<T> comparator
			)
		{
		final int indexes[]= equal_range(dataVector, first, last, select, comparator);
		if(indexes[0]==indexes[1]) return Collections.emptyIterator();
		return dataVector.subList(indexes[0], indexes[1]).iterator();
		}
	
	public static <T extends Comparable<T>> Stream<T> equal_range_stream(
			final List<T> dataVector,
			int first,
			final int last,
			final T select,
			final Comparator<T> comparator
			)
		{
		final Iterator<T> iter = equal_range_iterator(dataVector, first, last, select, comparator);
		return StreamSupport.stream(new IterableAdapter<T>(iter).spliterator(), false);
		}
	
}
