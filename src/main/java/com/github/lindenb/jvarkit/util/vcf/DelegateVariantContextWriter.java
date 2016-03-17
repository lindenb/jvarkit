package com.github.lindenb.jvarkit.util.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public class DelegateVariantContextWriter implements VariantContextWriter {
	private final VariantContextWriter _delegate;
	
	public DelegateVariantContextWriter(final VariantContextWriter delegate) {
		this._delegate = delegate;
	}
	
	public VariantContextWriter getDelegate() {
		return _delegate;
	}
	
	@Override
	public void add(final VariantContext ctx) {
		getDelegate().add(ctx);

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
	}

	@Override
	public void writeHeader(final VCFHeader header) {
		getDelegate().writeHeader(header);
	}

}
