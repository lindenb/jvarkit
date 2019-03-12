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
package com.github.lindenb.jvarkit.util.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** utility to extract allele frequencies from AC/AN or AF */
public class AFExtractorFactory {
	private static final Logger LOG = Logger.build(AFExtractorFactory.class).make();

/** an interface that can extract allele frequencies for each allele in a variant */
public static interface AFExtractor
	{
	/** validate the VCF header, return false on failure */
	public boolean validateHeader(final VCFHeader header);
	/** return the frequencies for each ALT allele. Item of the list may be null. List is never null */
	public List<Double> parse(final VariantContext ctx);
	/** return a now for this extractor */
	public String getName();
	}

public static final String OPT_DESC= 
	"How to extract the AlleleFrequencies from a variant. Multiple separated with comma or semicolon. e.g: \"AC/AN;exome_CEU_*;genome_NFE_AF;another_AC/another/AN\". Input is a set of AC/AN field pairs or/and AF field separated by semicolon. "
	+ "'x/y' means AC/AN fields. "
	+ "'*' will be replaced with AC and AN, hence, 'exome_CEU_*' will be interpreted as exome_CEU_AC/exome_CEU_AN. "
	+ "Other field will be interpreted as an AF field.";
/** par multiple AFExtractor, duplicate are removed. Syntax of user string : see OPT_DESC*/
public List<AFExtractor>  parseFieldExtractors(final String userStr) {
	final Set<AFExtractor> afset= new HashSet<>();
	if(!StringUtil.isBlank(userStr)) {
		for(final String s1:userStr.split("[;,]+")) {
			if(StringUtil.isBlank(s1)) continue;
			if(s1.contains("*") && s1.contains("/")) {
				throw new IllegalArgumentException("found both '/' and '*' in "+s1+" of "+userStr);
				}
			else if(s1.contains("*")) {
				afset.add(createAcAnFieldsExtractor(
						s1.replace("*","AC").trim(),
						s1.replace("*","AN").trim()
						));
				}
			else if(s1.contains("/")) {
				final int slash = s1.indexOf('/');
				afset.add(createAcAnFieldsExtractor(
						s1.substring(0,slash).trim(),
						s1.substring(slash+1).trim()
						));
				}
			else
				{
				afset.add(createAFFieldExtractor(s1));
				}
			}
		}
	return new ArrayList<>(afset);
	}


public AFExtractor createAFFieldExtractor(final String afAttr) {
	return new AFFieldExtractor(afAttr);
	}

private static class AFFieldExtractor implements AFExtractor
	{
	private String afAttr;
	AFFieldExtractor(final String afAttr) {
		this.afAttr = afAttr;
		if(StringUtil.isBlank(this.afAttr)) throw new IllegalArgumentException("AF is blank");
		}
	@Override
	public boolean validateHeader(final VCFHeader header) {
		final VCFInfoHeaderLine hdr = header.getInfoHeaderLine(this.afAttr);
		if(hdr==null) {
			final String msg = "VCF header doesn't contain the INFO field " +this.afAttr;
			LOG.warn(msg);
			return false;
			}
		if(hdr.getType()!=VCFHeaderLineType.Float) {
			final String msg = "INFO=<ID=" +this.afAttr+"> extected Type=Float but got "+hdr.getType();
			LOG.warn(msg);
			return false;
			}
		if(hdr.getCountType()!=VCFHeaderLineCount.A) {
			final String msg = "INFO=<ID=" +this.afAttr+"> extected Number=A but got "+hdr.getCountType();
			LOG.warn(msg);
			return false;
			}
		return true;
		}
	@Override
	public List<Double> parse(final VariantContext ctx) {
		if(!ctx.hasAttribute(this.afAttr)) return null;
		return ctx.getAttributeAsList(this.afAttr).stream().
				map(O->AFExtractorFactory.objectToFrequency(O)).
				collect(Collectors.toList());
		}
	
	@Override
	public int hashCode() {
		return this.afAttr.hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof AFFieldExtractor)) return false;
		return this.afAttr.equals(AFFieldExtractor.class.cast(obj).afAttr);
		}
	@Override
	public String getName() {
		return "INFO/"+this.afAttr;
		}
	@Override
	public String toString() {
		return "AF Extractor using INFO/"+this.afAttr;
		}
	}

public AFExtractor createAcAnFieldsExtractor(final String acAttr,final String anAttr) {
	return new ACANFieldsExtractor(acAttr,anAttr);
	}

private static final List<Double> double_list1 = Arrays.asList((Double)null);
private static final List<Double> double_list2 = Arrays.asList((Double)null,(Double)null);

private static class ACANFieldsExtractor implements AFExtractor
	{
	private final String acAttr;
	private final String anAttr;
	ACANFieldsExtractor(final String acAttr,final String anAttr) {
		this.acAttr = acAttr;
		this.anAttr = anAttr;
		if(StringUtil.isBlank(this.acAttr)) throw new IllegalArgumentException("AC is blank");
		if(StringUtil.isBlank(this.anAttr)) throw new IllegalArgumentException("AN is blank");
		if( (acAttr.contains("AN") && !acAttr.contains("AC")) && 
			(anAttr.contains("AC") && !anAttr.contains("AN"))) {
			LOG.warn("most probably the arguments have been swiched ac="+acAttr+";an="+anAttr);
			}
		if(this.acAttr.equals(this.anAttr)) throw new IllegalArgumentException(acAttr+"="+this.anAttr);
		}
	@Override
	public boolean validateHeader(final VCFHeader header) {
		VCFInfoHeaderLine hdr = header.getInfoHeaderLine(this.acAttr);
		if(hdr==null) {
			final String msg = "VCF header doesn't contain the INFO field " +this.acAttr;
			LOG.warn(msg);
			return false;
			}
		if(hdr.getType()!=VCFHeaderLineType.Integer) {
			final String msg = "INFO=<ID=" +this.acAttr+"> expected Type=Integer but got "+hdr.getType();
			LOG.warn(msg);
			return false;
			}
		if(hdr.getCountType()!=VCFHeaderLineCount.A) {
			final String msg = "INFO=<ID=" +this.acAttr+"> expected Number=A but got "+hdr.getCountType();
			LOG.warn(msg);
			return false;
			}
		hdr = header.getInfoHeaderLine(this.anAttr);
		if(hdr==null) {
			final String msg = "VCF header doesn't contain the INFO field " +this.anAttr;
			LOG.warn(msg);
			return false;
			}
		if(hdr.getType()!=VCFHeaderLineType.Integer) {
			final String msg = "INFO=<ID=" +this.anAttr+"> expected Type=Integer but got "+hdr.getType();
			LOG.warn(msg);
			return false;
			}
		
		if(!hdr.isFixedCount())
			{
			final String msg = "INFO=<ID=" +this.anAttr+"> expected fixed count "+hdr.getCountType();
			LOG.warn(msg);
			return false;
			}
		
		if(hdr.getCount()!=1) {
			final String msg = "INFO=<ID=" +this.anAttr+"> expected Number=1 but got "+hdr.getCount();
			LOG.warn(msg);
			return false;
			}
		return true;
		}
	
	
	
	private List<Double> emptyList(final VariantContext ctx) {
		final int n=ctx.getAlternateAlleles().size();
		switch(n)
			{
			case 1: return AFExtractorFactory.double_list1;
			case 2: return AFExtractorFactory.double_list2;
			default:
				final List<Double> L = new ArrayList<>(n);
				for(int i=0;i< n;++i) {
					L.add(null);
					}
				return L;
			}
		}
	
	@Override
	public List<Double> parse(final VariantContext ctx) {
		if(!ctx.hasAttribute(this.acAttr)) return emptyList(ctx);
		if(!ctx.hasAttribute(this.anAttr)) return emptyList(ctx);
		final double an;
		
		try {
		an = (double)ctx.getAttributeAsInt(this.anAttr, 0);
		} catch(final Throwable err) {
			LOG.error("bad value for INFO/"+this.anAttr+" at "+ctx.getContig()+"/"+ctx.getStart());
			return emptyList(ctx);
			}
		if(an>0) {
			return ctx.getAttributeAsList(this.acAttr).
					stream().
					map(O->AFExtractorFactory.objectToFrequency(O)).
					map(V->V==null?null:V/an).
					collect(Collectors.toList());
			}
		return emptyList(ctx);
		}
	
	@Override
	public int hashCode() {
		return this.acAttr.hashCode()*31+this.anAttr.hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof ACANFieldsExtractor)) return false;
		return this.acAttr.equals(ACANFieldsExtractor.class.cast(obj).acAttr) &&
				this.anAttr.equals(ACANFieldsExtractor.class.cast(obj).anAttr)
				;
		}
	@Override
	public String getName() {
		return "INFO/"+this.acAttr +"~INFO/"+this.anAttr;
		}
	
	@Override
	public String toString() {
		return "AF Extractor using " + getName();
		}
	}

private static Double objectToFrequency(final Object o) {
	if(o==null) return null;
	if(o instanceof Double) return Double.class.cast(o);
	if(o instanceof Float) return Float.class.cast(o).doubleValue();
	if(o instanceof Integer) return Integer.class.cast(o).doubleValue();
	final String s = o.toString();
	if(StringUtil.isBlank(s) || s.equals("null") || s.equals(".")) return null;
	return Double.parseDouble(s);
	}
}
