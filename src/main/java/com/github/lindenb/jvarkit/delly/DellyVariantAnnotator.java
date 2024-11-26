/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.delly;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
 * Why this class ? DELLY is missing INFO/SVLEN
 * @author lindenb
 *
 */
public class DellyVariantAnnotator implements VariantAnnotator {
	private boolean need_svlen;
	@Override
	public void fillHeader(final VCFHeader header) {
		need_svlen=header.getInfoHeaderLine("SVLEN")==null &&
				header.getInfoHeaderLine("CIPOS")!=null &&
				header.getInfoHeaderLine("CIEND")!=null &&
				header.getInfoHeaderLine("CT")!=null &&
				header.getInfoHeaderLine("IMPRECISE")!=null &&
				header.getInfoHeaderLine(VCFConstants.END_KEY)!=null &&
				header.getInfoHeaderLine("MAPQ")!=null &&
				header.getInfoHeaderLine("SVMETHOD")!=null;
		if(need_svlen) {
			header.addMetaDataLine(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer	,"SV Length"));
			}
		
		}
	
	public VariantContextBuilder fill(final VariantContextBuilder vcb,VariantContext src) throws IOException {
		if(need_svlen && src.hasAttribute(VCFConstants.END_KEY)) {
			final String svtype=src.getAttributeAsString("SVTYPE", "");
			if(svtype.equals("INV") || svtype.equals("DEL") || svtype.equals("DUP")) {
				final int svlen = src.getLengthOnReference() * (svtype.equals("DEL")?-1:1);
				vcb.attribute("SVLEN", svlen);
				}
			}
		return vcb;
		}
	
	public VariantContext fill(final VariantContext ctx) throws IOException {
		if(!need_svlen) return ctx;
		return fill(new VariantContextBuilder(ctx),ctx).make();
		}

	@Override
	public final List<VariantContext> annotate(VariantContext ctx) throws IOException {
		return Collections.singletonList(fill(ctx));
	}

	@Override
	public void close() {
	}

}
