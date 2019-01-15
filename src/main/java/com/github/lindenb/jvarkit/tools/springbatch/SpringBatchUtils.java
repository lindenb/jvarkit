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
package com.github.lindenb.jvarkit.tools.springbatch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.springframework.batch.item.ExecutionContext;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public class SpringBatchUtils {
public static final String VCF_HEADER_KEY="vcf-header";
public static VCFHeader getVcfHeader(final ExecutionContext executionContext) {
	if(!Objects.requireNonNull(executionContext, "executionContext NPE").containsKey(SpringBatchUtils.VCF_HEADER_KEY)) {
		throw new RuntimeException("key \""+SpringBatchUtils.VCF_HEADER_KEY+"\" is not defined in StepExecution");
		}
	final Object o =executionContext.get(SpringBatchUtils.VCF_HEADER_KEY);
	if(!(o instanceof VCFHeader)) {
		throw new RuntimeException("key \""+SpringBatchUtils.VCF_HEADER_KEY+"\" is not an instance of "+VCFHeader.class.getName()+" but "+o.getClass().getName());
		}
	return VCFHeader.class.cast(o);
	}

public static Map<String,Object> executionContextToMap(final ExecutionContext ctx) {
	if(ctx==null || ctx.isEmpty() ) return Collections.emptyMap();
	final Map<String,Object> hash = new HashMap<>(ctx.size());
	for(final Map.Entry<String, Object> key: ctx.entrySet()){
		hash.put(key.getKey(), key.getValue());
		}
	return hash;
	}

public static class VariantContextWriterBridge
	implements VariantContextWriter
	{
	private final List<VariantContext> variants = new ArrayList<VariantContext>();
	private VCFHeader header=null;
	
	@Override
	public void writeHeader(final VCFHeader header) {
		this.header = header;
		}
	public VCFHeader getHeader() {
		return this.header;
		}
	
	@Override
	public void add(final VariantContext vc) {
		this.variants.add(vc);
		}
	
	public void addAll(final List<VariantContext> vcs) {
		this.variants.addAll(vcs);
		}
	@Override
	public boolean checkError() {
		return false;
		}
	@Override
	public void close() {
		}
	public List<VariantContext> getAndClearVariants() {
		final List<VariantContext> L = new ArrayList<>(this.variants);
		this.variants.clear();
		return L;
		}
	}
}
