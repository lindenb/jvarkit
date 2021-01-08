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
package com.github.lindenb.jvarkit.samtools.util;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.Locatable;

/**
 * Basic implementation of Locatable for a single point position
 * @author lindenb
 *
 */
public class SimplePosition implements Locatable,Comparable<SimplePosition> {
	private final String contig;
	private final int pos;
	
	public SimplePosition(final String s) {
		if(s==null) throw new IllegalArgumentException("null argument");

		final int colon=s.indexOf(':');
		if(colon==-1 || colon+1==s.length())
			{
			throw new IllegalArgumentException("Bad 'chrom:pos' "+s);
			}
		
		this.contig=s.substring(0,colon).trim();
		if(StringUtils.isBlank(contig))
			{
			throw new IllegalArgumentException("Bad chrom:pos "+s);
			}
		this.pos=Integer.parseInt(s.substring(colon+1));
		}
	
	public SimplePosition(final String contig,final int pos) {
		this.contig = contig;
		if(this.contig==null) throw new IllegalArgumentException("contig is null");
		this.pos = pos;
		}

	@Override
	public String getContig() {
		return contig;
	}
	
	public int getPosition() {
		return this.pos;
	}
	
	@Override
	public int getStart() {
		return this.getPosition();
		}
	@Override
	public int getEnd() {
		return this.getPosition();
		}
	
	/** alias for getLengthOnReference */
	public final int length() {
		return 1;
	}
	
	@Override
	public int getLengthOnReference() {
		return 1;
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof SimplePosition)) return false;
 		final SimplePosition o = SimplePosition.class.cast(obj);
 		if(this.getPosition()!=o.getPosition()) return false;
		return this.getContig().equals(o.getContig());
		}
	
	@Override
	public int compareTo(final SimplePosition o) {
		int i= this.getContig().compareTo(o.getContig());
		if(i!=0) return i;
		i = Integer.compare(this.getPosition(),o.getPosition());
		return i;
		}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + contig.hashCode();
		result = prime * result + pos;
		return result;
	}
	
	@Override
	public String toString() {
		return this.contig+":"+this.pos;
	}
	
	public SimplePosition renameContig(final String ctg) {
		if(ctg.equals(this.getContig())) return this;
		return new SimplePosition(ctg,getPosition());
	}
	
	/** extends by dx bases */
	public Locatable extend(int dx) {
		if(dx==0) return this;
		if(dx>0) {
			return new SimpleInterval(getContig(),Math.max(1,getPosition()-dx),getPosition()+dx);
			}
		else
			{
			throw new IllegalArgumentException("negative extend");
			}
		}
	}
