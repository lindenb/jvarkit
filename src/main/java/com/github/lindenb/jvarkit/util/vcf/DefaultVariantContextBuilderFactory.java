package com.github.lindenb.jvarkit.util.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class DefaultVariantContextBuilderFactory
	implements VariantContextBuilderFactory {
public DefaultVariantContextBuilderFactory()
	{
	}
@Override
public VariantContextBuilder newVariantContextBuilder()
	{
	return new VariantContextBuilder();
	}
public VariantContextBuilder newVariantContextBuilder(final VariantContext ctx)
	{
	return new VariantContextBuilder(ctx);
	}
}
