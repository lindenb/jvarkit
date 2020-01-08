/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

public class HtsCollectors {

/** Filter Java Stream to 1 and only 1 element 
 * https://stackoverflow.com/questions/22694884/ 
 */
public static <T> Collector<T, ?, T> toSingleton() {
    return Collectors.collectingAndThen(
            Collectors.toList(),
            list -> {
            	switch(list.size())
            		{
            		case 0:  throw new IllegalStateException("expected at least one item but got none.");
            		case 1: return list.get(0);
            		default:  throw new IllegalStateException("expected only one item but got "+list.size());
            		}
            	}
    		);
		}

/** return item if there was one and only one item in the stream. 
 * If there is none or more than one, empty is returned.
 * Usage: want to take uniq item in set if set.size()==1
 *  */
public static <T> Collector<T, ?, Optional<T>> oneAndOnlyOne() {
    return Collectors.collectingAndThen(
            Collectors.toList(),
            list -> list.size()==1?Optional.of(list.get(0)):Optional.empty()
    		);
		}



/** convert stream of<QueryInterval> to an optimized array of QueryInterval */
public static Collector<QueryInterval, ?,QueryInterval[]>   optimizedQueryIntervals() {
	return Collectors.collectingAndThen(
		Collectors.toList(),
		list->{
			return  QueryInterval.optimizeIntervals(
					list.toArray(new QueryInterval[list.size()]));
			}
		);
	}

/*
public static <T extends Locatable> Collector<T, ?, Map<String,List<SimpleInterval>>>   mergeIntervals() {
	 return Collectors.collectingAndThen(
	            Collectors.groupingBy(R->R.getContig()),
	            mapIn -> {
	            	final Map<String,List<SimpleInterval>> listout= new HashMap<>(mapIn.size());
	            	for(final String contig: mapIn.keySet()) {
	            		final List<SimpleInterval> L = new ArrayList<>();
	            		for(Locatable item:mapIn.get(contig).stream().sorted(
			            		(A,B)-> Integer.compare(A.getStart(), B.getStart())
			            		).collect(Collectors.toList()))
	            			{
	            			if(!L.isEmpty() && L.get(L.size()-1).overlaps(item))
	            				{
	            				 L.set(L.size()-1,L.get(L.size()-1).merge(item));
	            				}
	            			else
	            				{
	            				L.add(new SimpleInterval(item));
	            				}
	            			}
	            		listout.put(contig, L);
	            		}
	            	return listout;
	            	}
	    		);
			}	
*/

public static <I extends Locatable> Collector<I,IntervalTreeMap<List<I>>,IntervalTreeMap<List<I>>>  
	toIntervalTreeMap()
	{
	return new DefaultCollector<>(
		()->new IntervalTreeMap<>(),
		(TMP,I)->{
			final Interval interval = new Interval(I.getContig(),I.getStart(),I.getEnd());;
			List<I> L = TMP.get(interval);
			if(L==null) {
				L = new ArrayList<>();
				TMP.put(interval,L);
				}
			L.add(I);
			},
		(left,right)->{
			for(final Interval interval:right.keySet())
				{
				List<I> L = left.get(interval);
				if(L==null) {
					L = new ArrayList<>();
					left.put(interval,L);
					}
				L.addAll(right.get(interval));
				}
			return left;
			}, 
		Collections.emptySet()
		);
	}


}
