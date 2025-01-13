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
package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;


import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.Lazy;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlEngine;
import org.apache.commons.jexl2.JexlException;
import org.apache.commons.jexl2.MapContext;
import org.apache.commons.jexl2.ObjectContext;


public class JexlBean implements Predicate<Object>, Function<Object,String> {
	private final static Lazy<JexlEngine> engine = new Lazy<>(() -> {
        final JexlEngine jexl = new JexlEngine();
        jexl.setSilent(false); // will throw errors now for selects that don't evaluate properly
        jexl.setLenient(false);
        jexl.setDebug(false);
        return jexl;
    });
	private Expression transform  = null;
	private String template = null;
	
	public JexlBean() {
		}
	public JexlBean(String template) {
		setTemplate(template);
	}
	public String getTemplate() {
		return template;
		}
	public void setTemplate(String template) {
		this.template = template;
		this.transform = null;
		}
	
	@Override
	public boolean test(Object t) {
		final Object o  = apply(t);
		if(o instanceof Boolean) return Boolean.class.cast(o).booleanValue();
		return false;
		}
	@Override
	public String apply(Object t) {
		return String.valueOf(eval(t));
		}
	
	
	protected Object eval(final Object o) {
		if(o==null || StringUtils.isBlank(getTemplate())) return "";
		try {
			if(this.transform==null) {
				this.transform = JexlBean.engine.get().createExpression(getTemplate());
				}
			final JexlContext context;
			if(o instanceof Map) {
				context = new MapContext(Map.class.cast(o));
				}
			else
				{
				context = new ObjectContext(JexlBean.engine.get(),o);
				}
			return transform.evaluate(context);
			}
		catch(JexlException err) {
			err.printStackTrace();
			return "";
			}
	}
}
