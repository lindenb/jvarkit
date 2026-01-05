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


import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Objects;
import java.util.Optional;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
 *  Base class for Analyzer
 *
 */
public abstract class AbstractAnalyzer implements Analyzer {
	protected static int ID_GENERATOR=0;
	protected Function<VCFHeader,Predicate<VariantContext>>  variantPredicateFactory = H->(V->true);
	protected Function<VCFHeader,Predicate<Genotype>>  genotypePredicateFactory = H->(G->true);
	private String analyzerName="";
	private String analyzerDesc="";
	private boolean enabled=true;
	protected VCFHeader vcfHeader;
	private Predicate<VariantContext>  variantPredicate =V->true;
	private Predicate<Genotype>  genotypePredicate = G->true;
	protected SamplePopulation sample2population = null;
	
	
	protected static class DoubleRounder implements DoubleUnaryOperator {
		final DecimalFormat decimalFormat;
		DoubleRounder(int n) {
			this("#."+StringUtils.repeat(n, '#'));
			}
		DoubleRounder(final String decimalFormatStr) {
			this.decimalFormat = new DecimalFormat(decimalFormatStr);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			}
		@Override
		public double applyAsDouble(double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}
		}
	
	protected static class BinRounder implements DoubleUnaryOperator {
		private final double[] bins;
		
		BinRounder(String commaStr) {
			this(Arrays.stream(CharSplitter.COMMA.split(commaStr)).mapToDouble(S->Double.parseDouble(S)).toArray());
			}
		BinRounder(double...values) {
			this.bins = Arrays.copyOf(values, values.length);
			Arrays.sort(this.bins);
			if(this.bins.length==0) throw new IllegalArgumentException();
			}
		@Override
		public double applyAsDouble(double v) {
			for(int i=0;i< this.bins.length;i++) {
				if(v <= this.bins[i]) return this.bins[i];
				}
			return this.bins[this.bins.length-1];
			}
		}
	
	
	protected AbstractAnalyzer() {
		this.analyzerName = getClass().getSimpleName().replace("Analyzer", "");
		this.analyzerDesc = analyzerName;
		}
	
	@Override
	public void init(VCFHeader vcfHeader, SamplePopulation sample2pop) {
		this.vcfHeader=vcfHeader;
		variantPredicate = Objects.requireNonNull(variantPredicateFactory.apply(vcfHeader));
		genotypePredicate = Objects.requireNonNull(genotypePredicateFactory.apply(vcfHeader));
		this.sample2population=sample2pop;
		}
	
	@Override
	public String getDescription() {
		return StringUtils.ifBlank(this.analyzerDesc,getName());
		}
	public void setDescription(String analyzerDesc) {
		this.analyzerDesc = analyzerDesc;
		}
	public void setName(String analyzerName) {
		this.analyzerName = analyzerName;
		}
	@Override
	public boolean isEnabled() {
		return this.enabled;
		}
	@Override
	public String getName() {
		return analyzerName;
		}
	
	@Override
	public void setVariantContextPredicateFactory(final Function<VCFHeader,Predicate<VariantContext>> factory) {
		this.variantPredicateFactory = factory;
		}
	
	public Predicate<VariantContext> getVariantContextPredicate() {
		return variantPredicate;
		}
	@Override
	public void setGenotypePredicateFactory(final Function<VCFHeader,Predicate<Genotype>> factory) {
		this.genotypePredicateFactory = factory;
		}
	public Predicate<Genotype> getGenotypePredicate() {
		return genotypePredicate;
		}
	
	protected boolean acceptVariant(final VariantContext ctx) {
		return getVariantContextPredicate().test(ctx);
		}
	protected boolean acceptGenotype(final Genotype ctx) {
		return getGenotypePredicate().test(ctx);
		}
	
	protected Optional<Genotype> getSingleton(VariantContext ctx) {
		Genotype single = null;
		for(Genotype g:ctx.getGenotypes()) {
			if(!acceptGenotype(g)) continue;
			if(!g.hasAltAllele()) continue;
			if(single!=null) return Optional.empty();
			single=g;
			}
		return Optional.ofNullable(single);
		}
	@Override
	public void disable() {
		this.enabled=false;
		}
}
