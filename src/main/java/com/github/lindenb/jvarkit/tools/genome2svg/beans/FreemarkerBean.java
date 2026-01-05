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
package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.io.StringWriter;
import java.util.function.Function;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.lang.StringUtils;

import freemarker.template.Configuration;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import htsjdk.samtools.util.Lazy;

public class FreemarkerBean implements Predicate<Object>, Function<Object,String> {
	private final static Lazy<Configuration> engine = new Lazy<>(() -> {
		final Configuration cfg = new  Configuration(Configuration.VERSION_2_3_32);
		return cfg;
    });

	
	private Template transform  = null;
	private String template = null;
	
	public FreemarkerBean() {
		}
	public FreemarkerBean(String template) {
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
		return Boolean.valueOf(apply(t));
		}
	@Override
	public String apply(Object t) {
		return eval(t);
		}
	
	
	public String eval(final Object o) {
		if(o==null || StringUtils.isBlank(getTemplate())) return "";
		try {
			if(this.transform==null) {
				this.transform = new Template("main",getTemplate(),FreemarkerBean.engine.get());
				}
			try(StringWriter sw = new StringWriter()) {
				this.transform.process(o, sw);
				return sw.toString();
				}
			}
		catch(IOException|TemplateException err) {
			err.printStackTrace();
			return "";
			}
	}
}
