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

import java.util.Objects;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public class DelegateVariantContextWriter implements VariantContextWriter {
	private final VariantContextWriter _delegate;
	private boolean is_closed=false;
	public DelegateVariantContextWriter(final VariantContextWriter delegate) {
		this._delegate = delegate;
	}
	
	public VariantContextWriter getDelegate() {
		return _delegate;
	}
	
	/** return true if 'close' was invoked */
	public boolean isClosed() {
		return this.is_closed;
	}
	
	@Override
	public void add(final VariantContext ctx) {
		if(isClosed()) return;
		Objects.requireNonNull(getDelegate()).add(ctx);
	}

	@Override
	public boolean checkError() {
		return getDelegate().checkError();
	}

	@Override
	public void close() {
		try { getDelegate().close(); } catch(final Throwable err) {
			/* https://github.com/samtools/htsjdk/issues/469 */
			}
		this.is_closed = true;
	}

	@Override
	public void setHeader(final VCFHeader header) {
		getDelegate().setHeader(header);
		}
	
	@Override
	public void writeHeader(final VCFHeader header) {
		getDelegate().writeHeader(header);
	}
	
	@Override
	public String toString() {
		return "instance-of(com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter)";
		}

}
