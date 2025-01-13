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
package com.github.lindenb.jvarkit.tools.vcfstats;

import com.github.lindenb.jvarkit.pedigree.SamplePopulation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SingletonAnalyzer extends AbstractFractionOfVariantsAnalyzer {
	public SingletonAnalyzer() {
		setName("Singletons");
		}
	@Override
	public void init(VCFHeader h, SamplePopulation sample2pop) {
		super.init(h, sample2pop);
		if(isEnabled() && h.getNGenotypeSamples()<=1) {
			disable();
			}
		}

	@Override
	public void visit(VariantContext ctx) {
		if(!acceptVariant(ctx)) return;
		super.n_variants++;
		final Genotype single = getSingleton(ctx).orElse(null);
		if(single!=null) super.counter.incr(single.getSampleName());
		}
	}
