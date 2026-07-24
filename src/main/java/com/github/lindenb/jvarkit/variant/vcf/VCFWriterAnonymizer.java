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
package com.github.lindenb.jvarkit.variant.vcf;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
/**
 * 
 * Anonymizer for VCF, well, it's a sample renamer too
 * 
 */
public class VCFWriterAnonymizer implements VariantContextWriter {
	private final VariantContextWriter delegate;
	private  final UnaryOperator<String> transform;
	private Map<String,String> remapper = null;
	public VCFWriterAnonymizer(final VariantContextWriter delegate,final UnaryOperator<String> transform) {
		this.delegate = delegate;
		this.transform = transform;
		}
	public VCFWriterAnonymizer(final VariantContextWriter delegate) {
		this(delegate,S->StringUtils.md5(S));
		}
	
	
	private VCFHeader renameHeader(final VCFHeader header) {
		if(!header.hasGenotypingData()) return header;
		if(this.remapper!=null) throw new IllegalStateException("remapped called twice");
		this.remapper=new HashMap<String, String>(header.getNGenotypeSamples());
		final Set<String> new_names= new HashSet<>(header.getNGenotypeSamples());
		for(String s: header.getSampleNamesInOrder()) {
			final String s2=  this.transform.apply(s);
			if(new_names.contains(s2)) throw new IllegalArgumentException("mapper produce the same name twice "+s2+" input was "+s);
			new_names.add(s2);
			this.remapper.put(s, s2);
			}
		return new VCFHeader(
				header.getMetaDataInInputOrder(),
				header.getSampleNamesInOrder().stream().map(S->this.remapper.get(S)).collect(Collectors.toList())
				);
		}
	
	@Override
	public void writeHeader(final VCFHeader header) {
		this.delegate.writeHeader(renameHeader(header));
	}

	@Override
	public void close() {
		this.delegate.close();		
	}

	@Override
	public boolean checkError() {
		return this.delegate.checkError();
	}

	@Override
	public void add(final VariantContext vc) {
		if(!vc.hasGenotypes()) {
			this.delegate.add(vc);
			}
		else
			{
			this.delegate.add(
				new VariantContextBuilder(vc)
					.genotypes(vc.getGenotypes().stream().map(G->new GenotypeBuilder(G).name(this.remapper.get(G.getSampleName())).make()).collect(Collectors.toList()))
					.make()
				);
			}
		}

	@Override
	public void setHeader(VCFHeader header) {
		this.delegate.setHeader(renameHeader(header));
	}

}
