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
package com.github.lindenb.jvarkit.jexl;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlException;
import org.apache.commons.jexl2.Script;

import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * A wrapper for org.apache.commons.jexl2.Expression / org.apache.commons.jexl2.Script
 */
public class BaseJexlHandler {
	public static final String OPT_WHAT_IS_JEXL="JEXL stands for Java EXpression Language.  See https://commons.apache.org/proper/commons-jexl/reference/syntax.html ";
	/** compiled expression */
	private final Expression jexlExpr;
	private final Script jexlScript;
	protected BaseJexlHandler(final String exprStr,boolean is_expression) {
		if(is_expression) {
			this.jexlExpr = VariantContextUtils.engine.get().createExpression(exprStr);
			this.jexlScript = null;
			}
		else
			{
			this.jexlExpr = null;
			this.jexlScript=VariantContextUtils.engine.get().createScript(exprStr);
			}
		}
	/** declared expression */
	public String getExpression() {
		return  isScript()?
				this.jexlScript.getText() :
				this.jexlExpr.getExpression();
		}

	
	protected boolean isExpression() {
		return this.jexlExpr!=null;
		}
	protected boolean isScript() {
		return !isExpression();
		}

	
	/** eval the context or throw on error */
	protected Object eval(final JexlContext ctx) {
		try {
			if(isScript())
				{
				return this.jexlScript.execute(ctx);
				}
			else
				{
				return this.jexlExpr.evaluate(ctx);
				}
			}
		catch(final JexlException err) {
			throw new RuntimeException("Cannot evaluate JEXL expression \""+this.getExpression()+"\".",err);
			}	
		}
	@Override
	public String toString() {
		return getClass().getName()+" "+this.getExpression();
		}
	
	}
