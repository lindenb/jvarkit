package com.github.lindenb.jvarkit.util.vcf;

import java.io.Closeable;
import java.util.Iterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/** net.sf.picard.vcf.VCFIterator deleted from from picard 1.100 */
public interface VcfIterator extends Iterator<VariantContext>,Closeable
	{
	public AbstractVCFCodec getCodec();
	public VCFHeader getHeader();
    public VariantContext peek();
    }
