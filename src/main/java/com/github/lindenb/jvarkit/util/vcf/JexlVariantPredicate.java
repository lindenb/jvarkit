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
package com.github.lindenb.jvarkit.util.vcf;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.PredicateVariantAnnotator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;

/** 
 * JEXL expression factory
 */
public class JexlVariantPredicate extends PredicateVariantAnnotator {
	private static final Logger LOG=Logger.of(JexlVariantPredicate.class);
	private static long ID_GENERATOR=System.currentTimeMillis();
	public static final String PARAMETER_DESCRIPTION = 
			"A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. " +
			"Empty string will accept all variants. " +
			"Expression returning a TRUE will accept the variant. "+
			"See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011 "
			;
	
	/*******************************************************************************************/
	private static final Predicate<VariantContext> ACCEPT_ALL=new Predicate<VariantContext>() {
		@Override
		public boolean test(VariantContext t) {
			return true;
		}
		@Override
		public String toString()
			{
			return "<empty string> (ACCEPT ALL)";
			}
	};
	/*******************************************************************************************/
	public static class Converter
	implements IStringConverter<Predicate<VariantContext> >
		{
		@Override
		public Predicate<VariantContext> convert(final String expr) {
			if(StringUtil.isBlank(expr)) return ACCEPT_ALL;
			return create(expr);
			}	
		}
	/*******************************************************************************************/
	
	public static Predicate<VariantContext> create(final String...exps) {
		return create(Arrays.asList(exps));
		}
	/*******************************************************************************************/

	public static Predicate<VariantContext> create(final List<String> exps) {
		final List<String> expressions = exps.stream().
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toList());
		if( expressions.isEmpty()) return ACCEPT_ALL;
		return new JexlVariantPredicate(exps);	
		}
	/*******************************************************************************************/
	private final List<JexlVCMatchExp> jexlVCMatchExps;
	
	public JexlVariantPredicate(final String expression) {
		this(Collections.singletonList(expression));
		}
	
	public JexlVariantPredicate(final List<String> expressions) {
		final List<String> dummyNames = expressions.stream().
				map(S->"JEXL"+(++ID_GENERATOR)).
				collect(Collectors.toList());
		try {
			this.jexlVCMatchExps = VariantContextUtils.initializeMatchExps(dummyNames, expressions);
			}
		catch(final Throwable err) {
			LOG.error(err);
			throw new RuntimeException("Cannot compile :"+String.join(",", expressions),err);
			}
		}
	
	@Override
	public boolean test(final VariantContext ctx) {
		try {
			return VariantContextUtils.match(ctx,this.jexlVCMatchExps).
				values().
				stream().
				anyMatch(B->B.booleanValue());
			}
		catch(final Throwable err) {
			throw new RuntimeException("JEXL Failed for variant "+ctx.getContig()+":"+ctx.getStart(),err);
			}
		}
	
	@Override
	public String getFilterDescription() {
		return "Filtered with JEXL expression(s) ( "+ PARAMETER_DESCRIPTION+" ) : "+ this.jexlVCMatchExps.stream().map(X->X.exp.getExpression()).collect(Collectors.joining(" "));
		}
	
	@Override
	public String toString() {
		return getClass().getName()+":"+
					this.jexlVCMatchExps.
					stream().
					map(x->x.toString()).
					collect(Collectors.joining(";"))
					;
		}
}
