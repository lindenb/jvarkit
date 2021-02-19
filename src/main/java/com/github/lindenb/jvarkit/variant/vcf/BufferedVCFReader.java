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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
 * A buffered VCFReader that stores the last query in memory
 *
 */
public class BufferedVCFReader implements VCFReader {
	public static final String OPT_BUFFER_DESC = "When we're looking for variant in a lare VCF file, load the variants in an interval of 'N' bases instead of doing a random access for each variant.";
	private final VCFReader delegate;
	private final int buffSizeInBp;
	private final List<VariantContext> buffer = new ArrayList<>();
	private Locatable lastInterval = null;
	private UnaryOperator<VariantContext> simplifier = V->V;
	private class MyIter extends AbstractCloseableIterator<VariantContext> {
		int i=0;
		final Locatable query;
		MyIter(Locatable query) {
			this.query = query;
			}
		@Override
		protected VariantContext advance() {
			while(i< buffer.size()) {
				final VariantContext ctx = buffer.get(i);
				i++;
				if(ctx.getStart()>this.query.getEnd()) {
					close();
					return null;
					}
				if(this.query.overlaps(ctx)) return ctx;
				}
			return null;
			}
		@Override
		public void close() {
			i=buffer.size();
			}
		}
	
	/** set a function to simplify (eg. remove genotypes) the variants. */
	public BufferedVCFReader setSimplifier(final UnaryOperator<VariantContext> simplifier) {
		this.simplifier = simplifier;
		return this;
		}
	
	/** 
	 * @param delegate the delegate {@link VCFReader}
	 * @param buffSizeInBp buffer size in bp
	 * */
	public BufferedVCFReader(final VCFReader delegate,int buffSizeInBp) {
		this.delegate = delegate;
		this.buffSizeInBp = buffSizeInBp;
		if(buffSizeInBp<1) throw new IllegalArgumentException("bad buffer size "+buffSizeInBp);
		}
	
	public VCFReader getDelegate() {
		return delegate;
		}
	
	/** close this and the delegate */
	@Override
	public void close() throws IOException {
		this.getDelegate().close();
		this.buffer.clear();
	}

	/* (non-Javadoc)
	 * @see htsjdk.variant.vcf.VCFReader#getHeader()
	 */
	@Override
	public VCFHeader getHeader() {
		return this.getDelegate().getHeader();
	}

	/** give a chance to remove genotypes, attributes, etc.. Variant will be ignored if returned value is null*/
	private VariantContext simplify(final VariantContext ctx) {
		return simplifier==null?ctx:simplifier.apply(ctx);
	}
	
	/* (non-Javadoc)
	 * @see htsjdk.variant.vcf.VCFReader#query(java.lang.String, int, int)
	 */
	@Override
	public CloseableIterator<VariantContext> query(final String chrom, int start, int end) {
		final Locatable query = new SimpleInterval(chrom,start,end);
		if(this.lastInterval==null || !this.lastInterval.contains(query)) {
			this.buffer.clear();
			this.lastInterval = new SimpleInterval(chrom, start, Math.max(end, start+this.buffSizeInBp));
			try(CloseableIterator<VariantContext> iter = this.getDelegate().query(this.lastInterval)) {
				while(iter.hasNext()) {
					final VariantContext ctx=simplify(iter.next());
					if(ctx==null) continue;
					this.buffer.add(ctx);
					}
				}
			}
		return new MyIter(query);
		}

	/* (non-Javadoc)
	 * @see htsjdk.variant.vcf.VCFReader#isQueryable()
	 */
	@Override
	public boolean isQueryable() {
		return this.getDelegate().isQueryable();
	}

	/* (non-Javadoc)
	 * @see htsjdk.variant.vcf.VCFReader#iterator()
	 */
	@Override
	public CloseableIterator<VariantContext> iterator() {
		return this.getDelegate().iterator();
	}

	@Override
	public String toString() {
		return this.getDelegate().toString();
		}
	
}
