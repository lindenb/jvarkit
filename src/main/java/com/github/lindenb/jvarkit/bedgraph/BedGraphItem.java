/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.bedgraph;


import java.util.Objects;

import com.github.lindenb.jvarkit.locatable.SimpleInterval;

import htsjdk.samtools.util.Locatable;

/**
 * an interval with a value
 */
public class BedGraphItem extends SimpleInterval {
	private final double value;
	public BedGraphItem(final Locatable loc,double value) {
		this(loc.getContig(),loc.getStart(),loc.getEnd(),value);
		}
	public BedGraphItem(final String contig,int start,int end,double value) {
		super(contig,start,end);
		this.value = value;
		}
	public double getValue() {
		return this.value;
		}
	
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj==null || !(obj instanceof BedGraphItem)) return false;
		final BedGraphItem o = (BedGraphItem) obj;
		if(this.getStart()!=o.getStart()) return false;
		if(this.getEnd()!=o.getEnd()) return false;
		if(!this.contigsMatch(o)) return false;
		return Double.doubleToLongBits(value) == Double.doubleToLongBits(o.value);
		}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + Objects.hash(value);
		return result;
		}
	
	@Override
	public String toString() {
		return getContig()+":"+getStart()+"-"+getEnd()+"="+getValue();
		}
	}
