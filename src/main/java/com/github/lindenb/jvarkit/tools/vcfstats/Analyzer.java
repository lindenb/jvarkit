/*
The MIT License (MIT)
Copyright (c) 2026 Pierre Lindenbaum

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


import java.io.IOException;
import java.nio.file.Path;
import java.util.function.Function;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.FunctionalMap;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
 *  Analyzer for VCFStats
 *
 */
public interface Analyzer {
	public void init(final VCFHeader h,final SamplePopulation sample2pop);
	public String getName();
	public String getDescription();
	public void setVariantContextPredicateFactory(final Function<VCFHeader,Predicate<VariantContext>> factory);
	public void setGenotypePredicateFactory(final Function<VCFHeader,Predicate<Genotype>> factory);
	public void visit(final VariantContext ctx);
	public boolean isEnabled();
	public void disable(); 
	public void writeReports(Path directory,final FunctionalMap<String, String> options) throws IOException;
}
