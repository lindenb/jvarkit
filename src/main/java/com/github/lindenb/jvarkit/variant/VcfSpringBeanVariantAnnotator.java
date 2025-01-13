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
package com.github.lindenb.jvarkit.variant;


import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.springframework.context.ApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VcfSpringBeanVariantAnnotator implements VariantAnnotator {
	public static final String DEFAULT_MAIN_BEAN_NAME="main";

private final List<String>  xmlConfigPaths = new ArrayList<>();
private String beanName= DEFAULT_MAIN_BEAN_NAME;
private final List<VariantAnnotator> annotators = new ArrayList<>();
private boolean initialized_flag = false;

public VcfSpringBeanVariantAnnotator() {
	// NO this(blablabala) !!
	}
public VcfSpringBeanVariantAnnotator(final Path xmlConfigPath,final String beanName) {
	this(Collections.singletonList(xmlConfigPath.toUri().toString()),beanName);
	}


public VcfSpringBeanVariantAnnotator(final List<String> xmlConfigs,final String beanName) {
	this.xmlConfigPaths.addAll(xmlConfigs);
	setBeanName(beanName);
	createAnnotators();
	}


public void setBeanName(String beanName) {
	this.beanName = beanName;
	}

public String getBeanName() {
	return beanName;
	}


private void createAnnotators() {
	if(this.initialized_flag) return;
	this.initialized_flag = true;
	if(this.xmlConfigPaths.isEmpty())
		{
		throw new IllegalArgumentException("no spring config file was provided");
		}
	try {
		
		final ApplicationContext springApplicationContext =  new FileSystemXmlApplicationContext(
				this.xmlConfigPaths.toArray(new String[this.xmlConfigPaths.size()])
				);
		
	
		
		if(!springApplicationContext.containsBean(getBeanName()))
			{
			throw new IllegalArgumentException("cannot get bean "+getBeanName()+" in "+ String.join(" ", this.xmlConfigPaths));
			}
		if(!springApplicationContext.containsBean(getBeanName())) {
			throw new IllegalArgumentException("bean "+getBeanName()+" is was not found in "+ String.join(" ", this.xmlConfigPaths));
			}
		final Object o = springApplicationContext.getBean(getBeanName());
		if( o == null) {
			throw new IllegalArgumentException("bean "+getBeanName()+" was found but it is null in "+ String.join(" ", this.xmlConfigPaths));
			}
		if( !(o instanceof List)) {
			throw new IllegalArgumentException("bean "+getBeanName()+" is not an instance of java.util.List in "+String.join(" ", this.xmlConfigPaths));
			}
		final List<?> list = (List<?>)o;
		for(int i=0;i< list.size();i++) {
			final Object o2 = list.get(i);
			if( o2 == null) {
				throw new IllegalArgumentException("element["+i+"]found but it is null in "+ String.join(" ", this.xmlConfigPaths));
				}
			if( !(o2 instanceof VariantAnnotator)) {
				throw new IllegalArgumentException("element["+i+"]found but it not an instance of "+VariantAnnotator.class.getName()+" in "+  String.join(" ", this.xmlConfigPaths));
				}
			this.annotators.add(VariantAnnotator.class.cast(o2));
			}
		}
	catch(final Throwable err) {
			throw new RuntimeException(err);
			}
	}

@Override
public void fillHeader(final VCFHeader header) {
	createAnnotators();
	for(VariantAnnotator ann:this.annotators) {
		ann.fillHeader(header);
		}
	}

private void recursive(int annotate_index, final VariantContext ctx, final List<VariantContext> filler) throws IOException {
	if(annotate_index==this.annotators.size()) {
		filler.add(ctx);
		}
	else
		{
		final VariantAnnotator ann = this.annotators.get(annotate_index);
		for(VariantContext ctx2: ann.annotate(ctx)) {
			recursive(annotate_index+1, ctx2, filler);
			}
		}
	}

@Override
public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
	final List<VariantContext> filler = new ArrayList<>();
	recursive(0,ctx,filler);
	return filler;
	}

@Override
public void close() {
	for(VariantAnnotator ann:this.annotators) {
		ann.close();
		}	
	}
}
