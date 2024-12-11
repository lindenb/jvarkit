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
import java.util.stream.Collector;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

public class LocatableUtils extends CoordMath {
	public static final Comparator<Locatable> DEFAULT_COMPARATOR = (A,B)->LocatableUtils.compareTo(A,B);
	
	public static Locatable parse(final String s,final SAMSequenceDictionary dictOrNull) {
		final Locatable loc0 = parse(s);
		if(dictOrNull==null) return loc0;
		String contig= loc0.getContig();
		SAMSequenceRecord ssr=dictOrNull.getSequence(contig);
		if(ssr==null) {
			if(contig.startsWith("chr")) {
				contig = contig.substring(3);
				}
			else
				{
				contig = "chr"+contig;
				}
			ssr=dictOrNull.getSequence(contig);
			if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(loc0.getContig(), dictOrNull);
			}
		else if(!ssr.getContig().equals(loc0.getContig()) /* alias */)
			{
			contig = ssr.getContig();
			}
		if(loc0.getStart()>ssr.getEnd()) throw new IllegalArgumentException("interval "+s+" is output chromosome "+ssr);
		return new SimpleInterval(contig,loc0.getStart(),loc0.getEnd());
		}
	
	public static Locatable parse(final String s) {
		final int col = s.lastIndexOf(":");
		if(col<=0) throw new IllegalArgumentException("no colon in "+s);
		final String contig = s.substring(0, col);
		final int hyphen = s.indexOf("-",col+1);
		int start;
		int end;
		if(hyphen<0)
			{
			start = Integer.parseInt(s.substring(col+1));
			end = start;
			}
		else
			{
			start= Integer.parseInt(s.substring(col+1, hyphen));
			end= Integer.parseInt(s.substring(hyphen+1));
			}
		if(start<0)  throw new IllegalArgumentException("start <0 in "+s);
		if(start>end)  throw new IllegalArgumentException("start>end in "+s);
		return new SimpleInterval(contig,start,end);
		}
	
	/** return the shared interval between two loc */
	public static Locatable sharedInterval(final Locatable loc1, final Locatable loc2) {
		if(!loc1.overlaps(loc2)) throw new IllegalArgumentException("no overlap between "+toString(loc1)+" and "+toString(loc2));
		return new SimpleInterval(loc1.getContig(),
				Math.max(loc1.getStart(), loc2.getStart()),
				Math.min(loc1.getEnd(), loc2.getEnd())
				);
		}	
	
	/** extend locatable to 'amount' fraction of interval length. Extends if amount > 1.0  */
	public static Locatable slopByFraction(final Locatable loc,double amount,final SAMSequenceDictionary dic) {
		return slopByBases(loc,(int)(loc.getLengthOnReference()*amount),dic);
		}
	
	/** extend locatable to 'amount' base pair. if amount==2, extends by 2 bp on 5' and 2 bp on 3'. shrink if amount <0  */
	public static Locatable slopByBases(final Locatable loc,int amount,final SAMSequenceDictionary dic) {
		final SAMSequenceRecord ssr = dic.getSequence(loc.getContig());
		if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), dic);
		if(loc.getStart() > ssr.getSequenceLength() || loc.getEnd()> ssr.getLengthOnReference()) throw new IllegalArgumentException(toString(loc)+" is out of bound of ref sequence "+toString(ssr));
		if(amount==0) {
			return loc;
			}
		else if(amount>0) {
			return new SimpleInterval(loc.getContig(),
					Math.max(1, loc.getStart()-amount),
					Math.min(loc.getEnd()+amount, ssr.getLengthOnReference())
					);
			}
		else
			{
			final int x1 = loc.getStart()+Math.abs(amount);
			final int x2 = loc.getEnd()-Math.abs(amount);
			if(x1==x2) return new  SimplePosition(loc.getContig(),x1);
			if(x1<=x2) return new  SimpleInterval(loc.getContig(),x1,x2);
			return new  SimplePosition(loc.getContig(),loc.getStart()+ loc.getLengthOnReference()/2);
			}
		}
	
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
	
	public static <T extends Locatable> Collector<T, ?, List<Locatable>> mergeIntervals() {
		return Collectors.collectingAndThen(Collectors.toList(), X->mergeIntervals(X));
		}
 	
	/** merge intervals on criteria. return a NEW list */
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
						Math.min(loc1.getStart(), loc2.getStart()),
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
	
	/** merge overlapping intervals if they are withing distance of 'x' . return a NEW list */
	public static List<Locatable> mergeIntervals(final List<? extends Locatable> L, final int distance) {
		return mergeIntervals(L,(A,B)->A.withinDistanceOf(B, distance));
		}
	
	/** merge overlapping intervals. return a NEW list */
	public static List<Locatable> mergeIntervals(final List<? extends Locatable> L) {
		return mergeIntervals(L,0);
		}
	
	
	/** merge two overlapping intervals. Do NOT check if they overlap, just check if they are on the same contig
	 * (allow merge is intervals are distant of 'x' bases */
	public static Locatable mergeIntervals(final Locatable loc1,final Locatable loc2) {
		if(!loc1.contigsMatch(loc2)) throw new IllegalArgumentException("intervals are not one the same contig "+toString(loc1)+" "+toString(loc2));
		return new SimpleInterval(
			loc1.getContig(),
			Math.min(loc1.getStart(), loc2.getStart()),
			Math.max(loc1.getEnd(),   loc2.getEnd())
			);
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
