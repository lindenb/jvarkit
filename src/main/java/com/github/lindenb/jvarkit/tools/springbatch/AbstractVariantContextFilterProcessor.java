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


*/
package com.github.lindenb.jvarkit.tools.springbatch;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public abstract class AbstractVariantContextFilterProcessor 
	extends AbstractVariantContextProcessor {	
public abstract String getFilterName();
public abstract Predicate<VariantContext> getPredicate();

@Override
public List<VariantContext> process(final List<VariantContext> variants) throws Exception {
	final Predicate<VariantContext> pred = this.getPredicate();
	final String filterName = this.getFilterName();
	if(pred == null || variants.isEmpty()) return variants;
	if(StringUtil.isBlank(filterName)) 
		{
		return variants.stream().
				filter(pred).
				collect(Collectors.toList());
		}
	
	return variants.stream().map(V->{
		if(pred.test(V)) {
			if(V.isFiltered()) return V;
			return new VariantContextBuilder(V).
					passFilters().
					make();
			}
		else
			{
			return new VariantContextBuilder(V).
					filter(filterName).
					make();
			}
		}).collect(Collectors.toList());
	}
}
