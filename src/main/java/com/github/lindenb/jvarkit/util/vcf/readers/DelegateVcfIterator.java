package com.github.lindenb.jvarkit.util.vcf.readers;

import java.io.IOException;

import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class DelegateVcfIterator implements VcfIterator{
	private final VcfIterator delegate;
	public DelegateVcfIterator(final VcfIterator delegate) {
		this.delegate=delegate;
		}
	
	protected VcfIterator getDelegate() {
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
	public void close() throws IOException {
		getDelegate().close();
		
	}

	@Override
	public AbstractVCFCodec getCodec() {
		return getDelegate().getCodec();
	}

	@Override
	public VCFHeader getHeader() {
		return getDelegate().getHeader();
	}

	@Override
	public VariantContext peek() {
		return getDelegate().peek();
	}

}
