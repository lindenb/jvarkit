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
package com.github.lindenb.jvarkit.variant.spring;

import java.util.function.Function;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
/**
 * class used for spring based application
 * keep variant using a jexl expression
 * @author lindenb
 *
 */
public class JEXLVariantPredicateFactory implements Function<VCFHeader,Predicate<VariantContext>>{
	private String expression;
	
	public JEXLVariantPredicateFactory() {
		this("true");
		}
	
	public JEXLVariantPredicateFactory(final String expression) {
		this.expression = expression;
		}
	
	public void setExpression(String expression) {
		this.expression = expression;
		}
	public String getExpression() {
		return expression;
		}
	@Override
	public Predicate<VariantContext> apply(final VCFHeader header) {
		return new JexlVariantPredicate(getExpression());		
		}
	}
