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
package com.github.lindenb.jvarkit.stream;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
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

/** convert a stream of locatable to a merged list of intervals */
public static  Collector<? super Locatable, ?,Stream<? extends Locatable>>   mergeIntervals() {
	return Collectors.collectingAndThen(
		Collectors.toCollection(ArrayList::new),
		list->{
			Collections.sort(list,(A,B)->{
				int i= A.getContig().compareTo(B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				i = Integer.compare(A.getEnd(), B.getEnd());
				return i;
				});
			int i=0;
			while(i +1 <list.size()) {
				final Locatable l1 = list.get(i);
				final Locatable l2 = list.get(i+1);
				if(l1.overlaps(l2)) {
					list.set(i,new SimpleInterval(
							l1.getContig(),
							Math.min(l1.getStart(), l2.getStart()),
							Math.max(l1.getEnd(), l2.getEnd())));
					list.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			return list.stream();
			}
		);
	}


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
