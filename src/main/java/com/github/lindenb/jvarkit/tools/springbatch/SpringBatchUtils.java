/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.springbatch;

import java.util.Objects;

import org.springframework.batch.item.ExecutionContext;

import htsjdk.variant.vcf.VCFHeader;

public class SpringBatchUtils {
public static final String VCF_HEADER_KEY="vcf-header";
public VCFHeader getVcfHeader(final ExecutionContext executionContext) {
	if(!Objects.requireNonNull(executionContext, "executionContext NPE").containsKey(SpringBatchUtils.VCF_HEADER_KEY)) {
		throw new RuntimeException("key \""+SpringBatchUtils.VCF_HEADER_KEY+"\" is not defined in StepExecution");
		}
	final Object o =executionContext.containsKey(SpringBatchUtils.VCF_HEADER_KEY);
	if(!(o instanceof VCFHeader)) {
		throw new RuntimeException("key \""+SpringBatchUtils.VCF_HEADER_KEY+"\" is not an instance of "+VCFHeader.class.getName()+" but "+o.getClass().getName());
		}
	return VCFHeader.class.cast(o);
	}
}
