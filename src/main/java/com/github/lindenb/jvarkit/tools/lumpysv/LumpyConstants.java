/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.lumpysv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Utilities for LUMPY-SV
 */
public class LumpyConstants {
	public static final Allele DEL = Allele.create("<DEL>", false);
	public static final Allele DUP = Allele.create("<DUP>", false);
	public static final Allele INV = Allele.create("<INV>", false);
	
	/** return true if the variant looks like a lumpy-sv header */
	public static boolean isLumpyHeader(final VCFHeader header) {
		final VCFFormatHeaderLine hL =header.getFormatHeaderLine("SU");
		return hL!=null;
		}
	/** return true if the variant looks like a lumpy-sv variant */
	public static boolean isLumpyVariant(final VariantContext ctx) {
		return	ctx.getStructuralVariantType()!=null &&
				ctx.getAlternateAlleles().size()==1 &&
				ctx.getAlternateAllele(0).isSymbolic() &&
				ctx.hasAttribute("SU")
				;
		}
}
