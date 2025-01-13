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
package com.github.lindenb.jvarkit.tools.vcfmulti2oneinfo;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** variant annotator converting multiple INFO par variant to multiple variant with one INFO */
public class MultiToOneInfoVariantAnnotator implements VariantAnnotator {
	private final String infoTag;
	public MultiToOneInfoVariantAnnotator(final String infoTag) {
		this.infoTag = infoTag;
		 if(StringUtil.isBlank(this.infoTag)) throw new IllegalArgumentException("empty INFO tag");
		}
	
	public String getTag() {
		return infoTag;
		}
	
	@Override
	public void fillHeader(final VCFHeader header) {
		final VCFInfoHeaderLine srcInfo = header.getInfoHeaderLine(getTag());
		if( srcInfo == null )
			{
			throw new IllegalArgumentException("Cannot find INFO FIELD 'INFO/"+ getTag()+"' in VCF header.");
			}
		switch( srcInfo.getCountType() )
			{
			case INTEGER:break;
			case UNBOUNDED:break;
			default: {
				throw new IllegalArgumentException("VCF header: INFO/CountType is not supported '"+ srcInfo.getCountType() +"'");
				}
			}
		switch( srcInfo.getType())
			{
			case Flag: 
				{
				throw new IllegalArgumentException("Type is not supported '"+ srcInfo.getType() +"'");
				}
			default:break;
			}
		}
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
		if(!ctx.hasAttribute(getTag())) {
			return Collections.singletonList(ctx);
			}
		
		final List<Object> L=ctx.getAttributeAsList(getTag());
		if( L.isEmpty() || L.size()==1)
			{
			return Collections.singletonList(ctx);
			}
		return L.stream().
			map(O->new VariantContextBuilder(ctx).attribute(getTag(),O).make()).
			collect(Collectors.toList());
		}
	
	@Override
	public void close() {
		}
	}
