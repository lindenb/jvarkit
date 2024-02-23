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
package com.github.lindenb.jvarkit.samtools.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiPredicate;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

public class LocatableUtils extends CoordMath {
	public static final Comparator<Locatable> DEFAULT_COMPARATOR = (A,B)->LocatableUtils.compareTo(A,B);
	
	/** convert locatable to simpleInterval or to simpleposition if length==1  */
	public static Locatable reduce(final Locatable A) {
		if(A.getStart()==A.getEnd()) {
			return new SimplePosition(A.getContig(), A.getStart());
			}
		else
			{
			return new SimpleInterval(A);
			}
		}
	/** return true of list is sorted using DEFAULT_COMPARATOR */
	public static boolean isSorted(List<? extends Locatable> L) {
		for(int i=0;i+1< L.size();++i) {
			if(DEFAULT_COMPARATOR.compare(L.get(i), L.get(i+1))>0) return false;
			}
		return true;
		}
	
	/** merge intervals on criteria */
	public static List<Locatable> mergeIntervals(final List<? extends Locatable> L0,final BiPredicate<Locatable, Locatable> predicate) {
		final List<Locatable> L=new ArrayList<>(L0);
		if(!isSorted(L)) Collections.sort(L,DEFAULT_COMPARATOR);
		int i=0;
		while(i +1 < L.size()) {
			final Locatable loc1 = L.get(i  );
			final Locatable loc2 = L.get(i+1);
			if(loc1.contigsMatch(loc2) && predicate.test(loc1,loc2)) {
				L.set(i,(Locatable)new SimpleInterval(
						loc1.getContig(),
						loc1.getStart(),
						Math.max(loc1.getEnd(), loc2.getEnd())
						));
				L.remove(i+1);
				}
			else
				{
				i++;
				}
			}
		return L;
		}
	
	/** merge overlapping intervals if they are withing distance of 'x' */
	public static List<Locatable> mergeIntervals(final List<? extends Locatable> L, final int distance) {
		return mergeIntervals(L,(A,B)->A.withinDistanceOf(B, distance));
		}
	
	/** merge overlapping intervals */
	public static List<Locatable> mergeIntervals(final List<? extends Locatable> L) {
		return mergeIntervals(L,0);
		}
	
	
	/** generic comparator on chrom/start/end */
	public static int compareTo(final Locatable A, Locatable B) {
		int i= A.getContig().compareTo(B.getContig());
		if(i!=0) return i;
		i = Integer.compare(A.getStart(),B.getStart());
		if(i!=0) return i;
		i = Integer.compare(A.getEnd(),B.getEnd());
		return i;
		}
	
	
	/** use chrom:start-end  */
	public static String toString(final Locatable loc) {
		return loc.getContig()+":"+ 
				String.valueOf(loc.getStart()) + "-"+ 
				String.valueOf(loc.getEnd())
				;
		}
	
	/** use StringUtils.niceInt to format start and end  */
	public static String toNiceString(final Locatable loc) {
		return loc.getContig()+":"+ 
				StringUtils.niceInt(loc.getStart()) + "-"+ 
				StringUtils.niceInt(loc.getEnd())
				;
		}
	/** use StringUtils.niceInt to format start and end  */
	public static String toBed3(final Locatable loc) {
		return String.join(
				"\t",
				loc.getContig(),
				String.valueOf(loc.getStart()-1) ,
				String.valueOf(loc.getEnd())
				);
		}
}
