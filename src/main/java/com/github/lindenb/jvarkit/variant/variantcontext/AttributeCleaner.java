/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant.variantcontext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** a class to remove fields from VariantContext similar to bcftools annotate */
public interface AttributeCleaner extends Function<VariantContext,VariantContext> {
public static final String OPT_DESC= "Variant Attribute cleaner. The syntax is the same as 'bcftools annotate'. e.g: 'INFO/AC,INFO/ANN'  Empty string does nothing.";

/** return true if vcf header should be kept */
public boolean isRetained(VCFHeaderLine hline);


public default Set<VCFHeaderLine> cleanMetaDataSet(final Set<VCFHeaderLine> set) {
	return set.stream().
				filter(H->isRetained(H)).
				collect(Collectors.toCollection(LinkedHashSet::new));
	}

public default VCFHeader cleanHeader(final VCFHeader header) {
	return new VCFHeader(cleanMetaDataSet(header.getMetaDataInInputOrder()),header.getSampleNamesInOrder());
	}

public VariantContextBuilder cleanVariant(VariantContext ctx,VariantContextBuilder vcb);

public default VariantContext cleanVariant(final VariantContext vc) {
	return cleanVariant(vc,new VariantContextBuilder(vc)).make();
}

@Override
public default VariantContext apply(final VariantContext ctx) {
	return cleanVariant(ctx);
}

/** compiles 'bctools annotate' like expression to a new AttributeCleaner */
public static AttributeCleaner compile(final String pattern) {
	if(pattern==null) return compile("");
	final AttributeCleanerImpl instance = new AttributeCleanerImpl(pattern);
	
	final String tokens[] = CharSplitter.COMMA.split(pattern);
	int i = 0;
	while(i< tokens.length) {
		String pat=tokens[i];
		if(StringUtils.isBlank(pat)) 
			{
			i++; continue;}
		if(pat.equals("ID")) {
			instance.remove_id = true;
			i++;
			continue;
			}
		if(pat.equals("FILTER")) {
			instance.remove_all_filters = true;
			i++;
			continue;
			}
		if(pat.equals("INFO")) {
			instance.remove_all_info = true;
			i++;
			continue;
			}
		if(pat.equals("FORMAT")) {
			instance.remove_all_format = true;
			i++;
			continue;
			}
		if(pat.equals("QUAL")) {
			instance.remove_qual = true;
			i++;
			continue;
			}
		boolean inverse = false;
		if(pat.startsWith("^"))
			{	
			inverse = true;
			pat = pat.substring(1);
			tokens[i]=pat;//replace tokens or we wont go in the while loop below
			if(pat.isEmpty()) throw new IllegalArgumentException("empty string '^'");
			}
		if(pat.startsWith("INFO/")) {
			final String prefix = "INFO/";
			instance.inverse_info = inverse;

			while(i< tokens.length && tokens[i].startsWith(prefix)) {
				final String field= tokens[i].substring(prefix.length());
				if(StringUtils.isBlank(field)) throw new IllegalArgumentException("bad format in "+pattern);
				instance.infos.add(field);
				i++;
				}
			i++;
			}
		else if(pat.startsWith("FORMAT/")) {
			instance.inverse_format= inverse;
			final String prefix = "FORMAT/";
			while(i< tokens.length && tokens[i].startsWith(prefix)) {
				final String field= tokens[i].substring(prefix.length());
				if(StringUtils.isBlank(field)) throw new IllegalArgumentException("bad format in "+pattern);
				if(!field.equals(VCFConstants.GENOTYPE_KEY)) {
					instance.formats.add(field);
					}
				i++;
				}
			i++;
			}
		else if(pat.startsWith("FILTER/")) {
			instance.inverse_filter= inverse;
			final String prefix = "FILTER/";
			while(i< tokens.length && tokens[i].startsWith(prefix)) {
				final String field= tokens[i].substring(prefix.length());
				if(StringUtils.isBlank(field)) throw new IllegalArgumentException("bad format in "+pattern);
				instance.filters.add(field);
				i++;
				}
			i++;
			}
		else
			{
			throw new IllegalArgumentException("bad format in "+pattern);
			}
		}
	
	return instance;
	}
	
static class AttributeCleanerImpl implements AttributeCleaner {
	final String pattern;
	boolean remove_all_info=false;
	boolean remove_all_filters = false;
	boolean remove_id  = false;
	boolean remove_qual  = false;
	boolean remove_all_format = false;
	boolean inverse_info= false;
	boolean inverse_format= false;
	boolean inverse_filter= false;
	Set<String> formats = new HashSet<>();
	Set<String> infos = new HashSet<>();
	Set<String> filters = new HashSet<>();
	
	AttributeCleanerImpl(final String pattern) {
		this.pattern = pattern;
		}
	
	@Override
	public int hashCode() {
		return this.pattern.hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof AttributeCleanerImpl)) return false;
		return this.pattern.equals(AttributeCleanerImpl.class.cast(obj).pattern);
		}
	
	@Override
	public boolean isRetained(final VCFHeaderLine hline) {
		if(hline instanceof VCFInfoHeaderLine) {
			if(remove_all_info) return false;
			final VCFInfoHeaderLine h = VCFInfoHeaderLine.class.cast(hline);
			boolean in_set = this.infos.contains(h.getID());
			return in_set == this.inverse_info;
			}
		else if(hline instanceof VCFFormatHeaderLine) {
			final VCFFormatHeaderLine h = VCFFormatHeaderLine.class.cast(hline);
			if(h.getID().equals(VCFConstants.GENOTYPE_KEY)) return true;
			if(remove_all_format) return false;
			boolean in_set = this.formats.contains(h.getID());
			return in_set == this.inverse_format;
			}
		else if(hline instanceof VCFFilterHeaderLine) {				
			if(remove_all_filters) return false;
			final VCFFilterHeaderLine h = VCFFilterHeaderLine.class.cast(hline);
			boolean in_set = this.filters.contains(h.getID());
			return in_set == this.inverse_filter;
			}
		return true;
		}
	
	@Override
	public VariantContextBuilder cleanVariant(final VariantContext ctx,final VariantContextBuilder vcb) {
		if(this.remove_qual) vcb.log10PError(VariantContext.NO_LOG10_PERROR);
		if(this.remove_id) vcb.noID();
		if(ctx.isFiltered() && (!this.filters.isEmpty() || this.remove_all_filters)) {
			if(this.remove_all_filters) {
				vcb.unfiltered();
			} else
				{
				final Set<String> set = new HashSet<>(ctx.getFilters());
				if(this.inverse_filter) {
					set.retainAll(this.filters);
				} else
					{
					set.removeAll(this.filters);
					}
				vcb.filters(set);
				}
		}
		
			
		if(this.remove_all_info)
			{
			vcb.attributes(Collections.emptyMap());
			}
		else if(!this.infos.isEmpty())
			{
			final Map<String,Object> atts = new HashMap<>(ctx.getAttributes());
			final Set<String> keys = new HashSet<>(atts.keySet());
			for(final String att: keys) {
				boolean remove = this.infos.contains(att);
				if(this.inverse_info) remove=!remove;
				if(remove) atts.remove(att);
				}
			vcb.attributes(atts);
			}
			
		if(ctx.hasGenotypes() && (this.remove_all_format || !this.formats.isEmpty()))
			{
			final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
			
			for(final Genotype gt:ctx.getGenotypes()) {
				final GenotypeBuilder gb=new GenotypeBuilder(gt);
				if(this.remove_all_format) {
					gb.noAD();
					gb.noPL();
					gb.noDP();
					gb.noGQ();
					gb.noAttributes();
					//gb.alleles(Collections.emptyList());
					}
				else
					{
					//internal attributes
					
					if(this.formats.contains(VCFConstants.GENOTYPE_PL_KEY) == !inverse_format) {
						gb.noPL();
						}
					if(this.formats.contains(VCFConstants.GENOTYPE_ALLELE_DEPTHS) == !inverse_format) {
						gb.noAD();
						}
					if(this.formats.contains(VCFConstants.DEPTH_KEY) == !inverse_format) {
						gb.noDP();
						}
					if(this.formats.contains(VCFConstants.GENOTYPE_QUALITY_KEY) == !inverse_format) {
						gb.noGQ();
						}
						
				    //extended attributes
					final Map<String,Object> atts = new HashMap<>(gt.getExtendedAttributes());
					final Set<String> keys = new HashSet<>(atts.keySet());
					for(final String key: keys) {
						boolean remove = this.formats.contains(key);
						if(inverse_format) remove=!remove;

						if(remove) atts.remove(key);
						}
					gb.noAttributes();
					gb.attributes(atts);
					}
				genotypes.add(gb.make());
			}
			
			vcb.genotypes(genotypes);
			}
		
		return vcb;
	}
	
	@Override
	public String toString() {
		return this.pattern;
		}
	}
}
