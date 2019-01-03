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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.util.AbstractMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
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
				ctx.hasAttribute("CIPOS") && 
				ctx.hasAttribute("CIEND")
				;
		}
	/*
	public static Interval getIntervalLeft(final VariantContext ctx)
		{
		if(!ctx.hasAttribute("CIPOS")) throw new IllegalArgumentException("No CIPOS in "+ctx);
		final List<Integer> ciposL= ctx.getAttributeAsIntList("CIPOS",0);
		if(ciposL.size()!=2) throw new IllegalArgumentException("len(CIPOS)!=2 in "+ctx);
		return new Interval(
			ctx.getContig(),
			ctx.getStart()+ciposL.get(0),
			ctx.getStart()+ciposL.get(1)
			);
		}
	
	public static Interval getIntervalRight(final VariantContext ctx)
		{
		if(!ctx.hasAttribute("CIEND")) throw new IllegalArgumentException("No CIEND in "+ctx);
		final List<Integer> ciendL= ctx.getAttributeAsIntList("CIEND",0);
		if(ciendL.size()!=2) throw new IllegalArgumentException("len(CIEND)!=2 in "+ctx);
		final String cL;
		final int pL;
		if(ctx.getStructuralVariantType()==StructuralVariantType.BND) {
			final  Map.Entry<String,Integer> entry = getBnDContigAndPos(ctx.getAlternateAllele(0).getDisplayString());
			cL = entry.getKey();
			pL = entry.getValue();
			} 
		else
			{
			cL = ctx.getContig();
			pL = ctx.getEnd();
			}
		
		return new Interval(
			cL,
			pL+ciendL.get(0),
			pL+ciendL.get(1)
			);
		}*/

	
	/** return Bnd contig and pos from a BnD ALT allele */
	public static Map.Entry<String,Integer> getBnDContigAndPos(final String alleleStr)
		{
		int left=-1,mid=-1,right=-1;
		int i=0;
		while(i< alleleStr.length())
			{
			char c= alleleStr.charAt(i);
			if(c=='[' ||c==']')
				{
				left=i;
				break;
				}
			i++;
			}
		if(left==-1) throw new IllegalArgumentException(alleleStr);

		while(i< alleleStr.length())
			{
			if(alleleStr.charAt(i)==':')
				{
				mid=i;
				break;
				}
			i++;
			}
		if(mid==-1) throw new IllegalArgumentException(alleleStr);
		while(i< alleleStr.length())
			{
			char c= alleleStr.charAt(i);
			if(c=='[' ||c==']')
				{
				right=i;
				break;
				}
			i++;
			}
		if(right==-1) throw new IllegalArgumentException(alleleStr);
		return new AbstractMap.SimpleEntry<String,Integer>(
				alleleStr.substring(left+1,mid),
				Integer.parseInt(alleleStr.substring(mid+1,right))
				);
		}
	
	public static String getBnDContig(final String alleleStr)
		{
		return getBnDContigAndPos(alleleStr).getKey();
		}
	public static int getBnDPos(final String alleleStr)
		{
		return getBnDContigAndPos(alleleStr).getValue();
		}
}
