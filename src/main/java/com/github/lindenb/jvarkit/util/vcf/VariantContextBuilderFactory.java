package com.github.lindenb.jvarkit.util.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public interface VariantContextBuilderFactory {
public VariantContextBuilder newVariantContextBuilder();
public VariantContextBuilder newVariantContextBuilder(final VariantContext ctx);
}
