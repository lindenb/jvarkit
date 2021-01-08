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
package com.github.lindenb.jvarkit.util.vcf.readers;


import htsjdk.variant.vcf.VCFIterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class DelegateVcfIterator implements VCFIterator{
	private final VCFIterator delegate;
	public DelegateVcfIterator(final VCFIterator delegate) {
		this.delegate=delegate;
		}
	
	protected VCFIterator getDelegate() {
		return delegate;
		}
	
	@Override
	public boolean hasNext() {
		return getDelegate().hasNext();
	}

	@Override
	public VariantContext next() {
		return getDelegate().next();
	}

	@Override
	public void close()  {
		getDelegate().close();
		
	}

	//@Override
	//public AbstractVCFCodec getCodec() { return getDelegate().getCodec();}

	@Override
	public VCFHeader getHeader() {
		return getDelegate().getHeader();
	}

	@Override
	public VariantContext peek() {
		return getDelegate().peek();
	}

}
