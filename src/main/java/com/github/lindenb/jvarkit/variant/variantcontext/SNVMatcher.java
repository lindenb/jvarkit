package com.github.lindenb.jvarkit.variant.variantcontext;
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
import java.util.Set;

import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

/**
 * Test if two variants overlap according to their ID, allele, etc...
 */
public enum SNVMatcher {
	id,
	overlap,
	chrom_pos,
	chrom_pos_ref,
	chrom_pos_ref_any_alt,
	chrom_pos_ref_all_alt
	;

	
	public static final String OPT_DESC = "How to test if two SNVs are the same. Spanning deletion will be ignored.";
	

	/** IMPORTANT: this function assumes that vc1 and vc2 are already known to overlap
	 * @param vc1 is the USER variant
	 * @param vc2 is the DATABASE variant
	 * @return true if the variants match using the matcher criteria
	 * 
	 * */
	public boolean test(final VariantContext vc1, final VariantContext vc2)
		{
		return test(this,vc1,vc2);
		}
	
	/** IMPORTANT: this function assumes that vc1 and vc2 are already known to overlap
	 * @param matcher the matcher
	 * @param vc1 is the USER variant
	 * @param vc2 is the DATABASE variant
	 * @return true if the variants match using the matcher criteria
	 * 
	 * */
	public static boolean test(
			final SNVMatcher matcher,
			final VariantContext vc1,
			final VariantContext vc2
			)
		{
		switch(matcher) {
			case overlap:
				return true;//assume we already know they overlap
			case id:
				{
				if(!vc1.hasID() || !vc2.hasID()) return false;
				final Set<String> id1s =CharSplitter.of(VCFConstants.ID_FIELD_SEPARATOR).splitAsSet((vc1.getID()));
				final Set<String> id2s =CharSplitter.of(VCFConstants.ID_FIELD_SEPARATOR).splitAsSet((vc2.getID()));
				for(String id: id1s) {
					if(id.equals(VCFConstants.EMPTY_ID_FIELD)) continue;
					if(id2s.contains(id)) return true;
					}
				return false;
				}
			case chrom_pos:
				return vc1.getStart() == vc2.getStart();
			case chrom_pos_ref:
				{
				if(vc1.getStart() != vc2.getStart()) return false;
				if(vc1.getNAlleles()<1 || vc2.getNAlleles()<1) return false;
				return vc1.getReference().equals(vc2.getReference());
				}
			case chrom_pos_ref_any_alt:
				{
				if(!test(SNVMatcher.chrom_pos_ref,vc1,vc2)) return false;
				for(Allele alt: vc1.getAlternateAlleles()) {
					if(alt.equals(Allele.SPAN_DEL)) continue;
					if(vc2.hasAllele(alt)) return true;
					}
				return false;
				}
			case chrom_pos_ref_all_alt:
				{
				if(!test(SNVMatcher.chrom_pos_ref,vc1,vc2)) return false;
				if(vc1.getNAlleles()<2) return false;
				if(vc2.getNAlleles()<2) return false;
				boolean ok=false;
				for(Allele alt: vc1.getAlternateAlleles()) {
					if(alt.equals(Allele.SPAN_DEL)) continue;
					if(!vc2.hasAllele(alt)) return false;
					ok=true;
					}
				return ok;
				}
			default: throw new IllegalStateException(matcher.name());
			}
		}
}
