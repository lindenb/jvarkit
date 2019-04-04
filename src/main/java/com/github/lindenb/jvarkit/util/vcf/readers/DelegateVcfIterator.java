package com.github.lindenb.jvarkit.util.vcf.readers;

import java.io.IOException;

import htsjdk.variant.vcf.VCFIterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
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
