/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.locatable;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

/**
 * Basic implementation of Locatable
 * @author lindenb
 *
 */
public class SimpleInterval extends AbstractLocatable implements Comparable<SimpleInterval> {
	private final String contig;
	private final int start;
	private final int end;
	
	public SimpleInterval(final String s) {
		if(s==null) throw new IllegalArgumentException("null argument");
		final int colon =s.indexOf(":");
		if(colon==-1) throw new IllegalArgumentException("no colon in "+s);
		this.contig = s.substring(0,colon);
		final String s2=s.substring(colon+1).replace(",", "");
		final int hyphen = s2.indexOf('-');
		if(hyphen==-1)
			{
			final int plus= s2.indexOf('+');
			if(plus==-1) {
				this.start = Integer.parseInt(s.substring(colon+1));
				this.end=this.start;
				}
			else
				{
				final int pivot = Integer.parseInt(s2.substring(0,plus));
				final int extend = Integer.parseInt(s2.substring(plus+1));
				this.start = Math.max(1,pivot-extend);
				this.end = pivot+extend;
				}
			}
		else
			{
			this.start = Integer.parseInt(s2.substring(0,hyphen));
			this.end = Integer.parseInt(s2.substring(hyphen+1));
			}
		}
	public SimpleInterval(final Locatable loc) {
		this(loc.getContig(),loc.getStart(),loc.getEnd());
	}
	
	/** no this is ambigous 
	public SimpleInterval(final BedInterval loc) {
		this(loc.getContig(),loc.getBedStart()+1,loc.getBedEnd());
	} */
	
	public SimpleInterval(final String contig,final int start,final int end) {
		this.contig = contig;
		if(this.contig==null) throw new IllegalArgumentException("contig is null");
		this.start=start;
		this.end = end;
	}

	@Override
	public String getContig() {
		return contig;
	}
	@Override
	public int getStart() {
		return start;
		}
	@Override
	public int getEnd() {
		return end;
		}

	
	
	/** alias for getLengthOnReference */
	public final int length() {
		return this.getLengthOnReference();
	}
	
	@Override
	public boolean equals(final Object obj) {
		return super.equals(obj);
		}
	
	@Override
	public int compareTo(final SimpleInterval o) {
		return LocatableUtils.compareTo(this, o);
		}
	
	
	
	public SimpleInterval renameContig(final String ctg) {
		if(ctg.equals(this.getContig())) return this;
		return new SimpleInterval(ctg,getStart(),getEnd());
	}
	
	/** get the merged interval <code>SimpleInterval(chrom,min(this.getStart(),other.getStart()),max(this.getStart(),other.getStart()))</code> */
	public SimpleInterval merge(final Locatable other) {
		if(!overlaps(other)) throw new IllegalArgumentException("assertion failed "+this.toNiceString() + " doesn't overlap with "+other);
		if(this.contains(other)) return this;
		return new SimpleInterval(
				getContig(),
				Math.min(getStart(),other.getStart()),
				Math.max(getEnd(),other.getEnd())
				);
		}
	
	public SimpleInterval extend(int dx) {
		if(dx==0) return this;
		if(dx>0) {
			return new SimpleInterval(getContig(),Math.max(1,getStart()-dx),getEnd()+dx);
			}
		else
			{
			throw new IllegalArgumentException("negative extend");
			}
		}
	
	public static List<SimpleInterval> mergeCollection(final Collection<? extends Locatable> locs) {
		final List<SimpleInterval> r= locs.stream().
				map(R->new SimpleInterval(R)).
				sorted().
				collect(Collectors.toCollection(ArrayList::new))
				;
		int i=0;
		while(i+1 < r.size()) {
			final SimpleInterval r1 = r.get(i);
			final SimpleInterval r2 = r.get(i+1);
			if(r1.overlaps(r2)) {
				r.set(i, r1.merge(r2));
				r.remove(i+1);
				}
			else
				{
				i++;
				}
			}
		return r;
		}
	/** get center position */
	public SimplePosition getCenter() {
		return new SimplePosition(getContig(), getStart()+this.getLengthOnReference()/2);
		}
	
	/** return intersection length. 0 if not same contig */
	 public int getIntersectionLength(final Locatable other) {
     if (this.overlaps(other)) {
         return CoordMath.getOverlap(this.getStart(), this.getEnd(), other.getStart(), other.getEnd());
     	}
     return 0;
	 }
	}
