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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** variant annotator converting multiple ALT par variant to multiple variant with one ALT */
public class MultiToOneAllele implements Function<VariantContext,List<VariantContext>> {
	private final Set<VCFInfoHeaderLine> infoHeaderToChange;
	private final Map<String,VCFFormatHeaderLine> formatHeaderToChange;
	private boolean keepOnlyMostFrequentAltAllele = false;

	public MultiToOneAllele(final VCFHeader header) {
		this.infoHeaderToChange = header.getInfoHeaderLines().stream().
			filter(H->H.getCountType().equals(VCFHeaderLineCount.A) || H.getCountType().equals(VCFHeaderLineCount.R)).
			collect(Collectors.toSet());
		this.formatHeaderToChange = header.getFormatHeaderLines().stream().
				filter(H->H.getCountType().equals(VCFHeaderLineCount.A) || H.getCountType().equals(VCFHeaderLineCount.R)).
				collect(Collectors.toMap(H->H.getID(),H->H));
		}
	
	public MultiToOneAllele setKeepOnlyMostFrequentAltAllele(boolean keepOnlyMostFrequentAltAllele) {
		this.keepOnlyMostFrequentAltAllele = keepOnlyMostFrequentAltAllele;
		return this;
		}
	
	
	
	@Override
	public List<VariantContext> apply(final VariantContext ctx)  {
		final List<Allele> alternateAlleles = ctx.getAlternateAlleles();
		if(alternateAlleles.isEmpty()) return Collections.emptyList();
		if(alternateAlleles.size()==1) return Collections.singletonList(ctx);
		
		final List<VariantContext> variants = new ArrayList<>(alternateAlleles.size());
		final Allele mostFrequent = this.keepOnlyMostFrequentAltAllele?ctx.getAltAlleleWithHighestAlleleCount():null;
		
		
		for(int alternateIndex=0;alternateIndex< alternateAlleles.size();++alternateIndex)
			{
			final Allele the_alt = alternateAlleles.get(alternateIndex);
			if(mostFrequent!=null && !the_alt.equals(mostFrequent)) continue;
			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.alleles(Arrays.asList(ctx.getReference(),the_alt));				
			
			for(VCFInfoHeaderLine header: this.infoHeaderToChange) {
				if(!ctx.hasAttribute(header.getID())) continue;
				final List<Object> list = ctx.getAttributeAsList(header.getID());
				
				if(header.getCountType().equals(VCFHeaderLineCount.R) && alternateIndex+1 < list.size())
					{
					vcb.attribute(header.getID(), Arrays.asList(
							list.get(0)/* REF */,
							list.get(alternateIndex+1))
							);	
					}
				else if(header.getCountType().equals(VCFHeaderLineCount.A) && alternateIndex < list.size())
					{	
					vcb.attribute(header.getID(), list.get(alternateIndex));	
					}
				else
					{
					vcb.rmAttribute(header.getID());
					}
				}
			
			if(ctx.hasGenotypes()) {
				final List<Genotype> genotypes=new ArrayList<>(ctx.getNSamples());
				for(int i=0;i< ctx.getNSamples(); i++) {
					final Genotype g= ctx.getGenotype(i);
					final GenotypeBuilder gb =new GenotypeBuilder(g);
					
					final List<Allele> gtAlleles = new ArrayList<>(g.getAlleles());
					for(int j=0;j< gtAlleles.size();j++) {
						Allele a = gtAlleles.get(j);
						if(a.isNoCall() || a.isReference() || a.equals(the_alt)) continue;
						gtAlleles.set(j, ctx.getReference());
						}
					gb.alleles(gtAlleles);
					gb.noAttributes();
					
					if(g.hasAD()) gb.AD(reduceArray(g.getAD(),alternateIndex+1));
					if(g.hasPL()) gb.PL(reduceArray(g.getPL(),alternateIndex+1));
					
					for(final String key:g.getExtendedAttributes().keySet()) {
						final VCFFormatHeaderLine header = this.formatHeaderToChange.getOrDefault(key, null);
						if(header==null) {
							gb.attribute(key, g.getExtendedAttribute(key));
							continue;
							}
						
						
						final Object oL = g.getExtendedAttribute(key);
						if(oL==null || !(oL instanceof List)) throw new IllegalArgumentException("not a List but" +oL);
						final List<?> list = (List<?>)oL;
						
						if(header.getCountType().equals(VCFHeaderLineCount.R) && alternateIndex+1 < list.size())
							{
							gb.attribute(header.getID(), Arrays.asList(
									list.get(0)/* REF */,
									list.get(alternateIndex+1))
									);	
							}
						else if(header.getCountType().equals(VCFHeaderLineCount.A) && alternateIndex < list.size())
							{	
							gb.attribute(header.getID(), list.get(alternateIndex));	
							}
						else
							{
							throw new IllegalArgumentException(key+" for "+g.getSampleName()+" in "+ctx.getContig()+":"+ctx.getStart());
							}
						}
					
					genotypes.add(gb.make());
					}
				vcb.genotypes(genotypes);
				}
			
			variants.add(vcb.make());
			}
		return variants;
		}
	
	private int[] reduceArray(int[] array,int keepIndex) {
		final int[] array2=new int[array.length-1];
		int i=0;
		int j=0;
		for(i=0;i< array.length;i++) {
			if(i==0 || i==keepIndex) {
				array2[j] = array[i];
				j++;
				}
			}
		return array2;
		}

	}
