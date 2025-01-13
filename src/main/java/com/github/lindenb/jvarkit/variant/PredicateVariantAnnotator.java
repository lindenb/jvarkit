/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

/** abstract <code>Predicate<VariantContext></code> that acts like a Predicate<VariantContext> */
public abstract class PredicateVariantAnnotator implements VariantAnnotator,Predicate<VariantContext> {
private boolean soft_filtering = true;
private VCFFilterHeaderLine filterHeaderLine = null;
private String filterName;

public PredicateVariantAnnotator() {
	this.filterName  = getClass().getSimpleName();
	}

public void setSoftFiltering(boolean b) {
	this.soft_filtering = b;
	}
public boolean isSoftFiltering() {
	return soft_filtering;
	}
/** filter name if soft filtering */
public String getFilterName() {
	return this.filterName;
	}

public void setFilterName(String filterName) {
	this.filterName = filterName;
	}

/** filter description if soft filtering */
public String getFilterDescription() {
	return "Filtered with "+getFilterName();
	}
@Override
public abstract boolean test(VariantContext ctx);
@Override
public void fillHeader(final VCFHeader header) {
	if(isSoftFiltering()) {
		if(StringUtils.isBlank(getFilterName())) throw new IllegalArgumentException("FILTER name is blank");
		this.filterHeaderLine = new VCFFilterHeaderLine(
			getFilterName(),
			getFilterDescription()
			);
		header.addMetaDataLine(this.filterHeaderLine);
	} else {
			this.filterHeaderLine = null;
		}
	}

@Override
public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
	final boolean accept= test(ctx);
	if(this.filterHeaderLine!=null) /* soft filtering */ {
		if(accept) {
			if(ctx.isFiltered()) return Collections.singletonList(ctx);
			return Collections.singletonList(new VariantContextBuilder(ctx).passFilters().make());
			}
		else
			{
			return Collections.singletonList(new VariantContextBuilder(ctx).filter(this.filterHeaderLine.getID()).make());
			}
		}
	else
		{
		return accept?Collections.singletonList(ctx):Collections.emptyList();
		}
	}

@Override
public void close() {
	}
}
