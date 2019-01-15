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
package com.github.lindenb.jvarkit.util.vcf;

import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;

/** 
 * JEXL expression factory
 */
public class JexlVariantPredicate implements Predicate<VariantContext> {
	private static final Logger LOG=Logger.build(JexlVariantPredicate.class).make();
	private static long ID_GENERATOR=System.currentTimeMillis();
	public static final String PARAMETER_DESCRIPTION = 
			"A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. " +
			"Empty string will accept all variants. " +
			"Expression returning a TRUE will accept the variant. "+
			"See https://gatkforums.broadinstitute.org/gatk/discussion/1255 "
			;
	
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
	
	public static class Converter
	implements IStringConverter<Predicate<VariantContext> >
		{
		@Override
		public Predicate<VariantContext> convert(final String expr) {
			if(StringUtil.isBlank(expr)) return ACCEPT_ALL;
			return create(expr);
			}	
		}

	
	public static Predicate<VariantContext> create(final String...exps) {
		return create(Arrays.asList(exps));
		}
	public static Predicate<VariantContext> create(final List<String> exps) {
		final List<String> expressions = exps.stream().
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toList());
		if( expressions.isEmpty()) return ACCEPT_ALL;
		
		final List<String> dummyNames = expressions.stream().
				map(S->"JEXL"+(++ID_GENERATOR)).
				collect(Collectors.toList());
		try {
			return new JexlVariantPredicate(VariantContextUtils.initializeMatchExps(dummyNames, expressions));
			}
		catch(final Throwable err) {
			LOG.error(err);
			throw new RuntimeException("Cannot compile :"+String.join(",", expressions),err);
			}
	}
	
	private final List<JexlVCMatchExp> jexlVCMatchExps;
	
	private JexlVariantPredicate(final List<JexlVCMatchExp> jexlVCMatchExps) {
		this.jexlVCMatchExps = jexlVCMatchExps;
		if(jexlVCMatchExps==null) throw new RuntimeException("jexlVCMatchExps is null");
		}
	@Override
	public boolean test(final VariantContext ctx) {
		return VariantContextUtils.match(ctx,this.jexlVCMatchExps).
			values().
			stream().
			anyMatch(B->B.booleanValue());
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
