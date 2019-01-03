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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.vcf;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * 
 * Class wrapping chrom/pos/ref, to be used as a key in a Map<ContigPosRef,something>
 *
 */
public class ContigPosRef implements htsjdk.samtools.util.Locatable,Comparable<ContigPosRef>{
	private final String contig;
	private final int pos;
	private final Allele ref;
	public ContigPosRef(final String contig,final int pos,final Allele ref)
		{
		this.contig = contig;
		this.pos = pos;
		this.ref = ref;
		}
	
	public ContigPosRef(final VariantContext ctx)
		{
		this(ctx.getContig(),ctx.getStart(),ctx.getReference());
		}
	public boolean isSameAs(final VariantContext ctx) {
		if( ctx == null ) return false;
		return new ContigPosRef(ctx).equals(this);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + contig.hashCode();
		result = prime * result + pos;
		result = prime * result + ref.hashCode();
		return result;
	}

	@Override
	public int compareTo(final ContigPosRef o) {
		int i = this.contig.compareTo(o.contig);
		if( i!=0 ) return i; 
		i = this.pos - o.pos;
		if( i!=0 ) return i; 
		return this.ref.compareTo(o.ref);
		}
	
	@Override
	public boolean equals(final Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final ContigPosRef other = (ContigPosRef) obj;
		return  this.pos==other.pos &&
				this.contig.equals(other.contig) &&
				this.ref.equals(other.ref);
	}

	@Override
	public String getContig() {
		return contig;
		}
	public int getPos() {
		return pos;
		}
	public Allele getReference() {
		return ref;
		}
	/** same as getPos  */
	@Override
	public int getStart() {
		return getPos();
		}
	/** getStart()+ getReference().length()-1 */
	@Override
	public int getEnd() {
		return getStart()+getReference().length()-1;
		}
	@Override
	public String toString() {
		return getContig()+":"+getPos()+"("+getReference().getDisplayString()+")";
		}
	}
