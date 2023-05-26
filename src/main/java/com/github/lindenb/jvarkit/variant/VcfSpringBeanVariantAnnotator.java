/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.List;

import org.springframework.context.ApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VcfSpringBeanVariantAnnotator implements VariantAnnotator {
private final List<String> springCongigFiles = new ArrayList<>();
private final List<VariantAnnotator> annotators = new ArrayList<>();
public VcfSpringBeanVariantAnnotator() {
	}


public void setConfigFiles(final List<String> springCongigFiles) {
	if(springCongigFiles==null) throw new IllegalArgumentException("springCongigFiles==nunll");
	this.springCongigFiles.clear();
	this.springCongigFiles.addAll(springCongigFiles);
	}

private void createAnnotators() {
	String mainBeanName = "main";
	if(this.springCongigFiles.isEmpty())
		{
		throw new IllegalArgumentException("no spring config file was provided");
		}
	final ApplicationContext springApplicationContext = 
			new FileSystemXmlApplicationContext(
					this.springCongigFiles.toArray(new String[this.springCongigFiles.size()])
					);
	
	if(!springApplicationContext.containsBean(mainBeanName))
		{
		throw new IllegalArgumentException("cannot get bean "+mainBeanName+" in "+
				String.join(" ",springCongigFiles)
				); 
		}
	final Object o = springApplicationContext.getBean(mainBeanName);
	if( o == null || !(o instanceof List))
		{
		throw new IllegalArgumentException("bean "+mainBeanName+" is not a  VcfChain but "+
				(o==null?"null":o.getClass().getName())
				); 
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
