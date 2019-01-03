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
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/** recalculate variant context attribute AF, AC, AN */
public class VariantAttributesRecalculator {
@Parameter(names={"--disable-vc-attribute-recalc"},description="When genotypes are removed/changed, Dd not recalculate variant attributes like DP, AF, AC, AN...")
private boolean recalculation_is_disabled = false;
@Parameter(names={"--vc-attribute-recalc-ignore-filtered"},description="When recalculating variant attributes like DP AF, AC, AN, ignore FILTERed **Genotypes**")
private boolean recalculation_ignore_gt_filtered = false;
@Parameter(names={"--vc-attribute-recalc-ignore-missing"},description="Ignore missing VCF headers (DP, AF, AC, AN). Default behavior: adding VCF header if they're missing")
private boolean do_not_add_missing_headers = false;



private boolean header_was_set=false;
private boolean DP_flag=false;
private boolean AC_flag=false;
private boolean AN_flag=false;
private boolean AF_flag=false;

/** add missing field to the header */
public VCFHeader setHeader(final VCFHeader header) {
	
	this.header_was_set = true;
	this.DP_flag = header.getInfoHeaderLine(VCFConstants.DEPTH_KEY)!=null;
	this.AC_flag = header.getInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY)!=null;
	this.AN_flag = header.getInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY)!=null;
	this.AF_flag = header.getInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY)!=null;
	if(!do_not_add_missing_headers) {
		final Set<VCFHeaderLine> headerLines = new HashSet<>();
		if(!this.DP_flag) {
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DEPTH_KEY);
			//VCFStandardHeaderLines.addStandardFormatLines(headerLines, true, VCFConstants.DEPTH_KEY);
		}
		if(!this.AC_flag) VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.ALLELE_COUNT_KEY);
		if(!this.AN_flag) VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.ALLELE_NUMBER_KEY);
		if(!this.AF_flag) VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.ALLELE_FREQUENCY_KEY);
		this.AC_flag = true;
		this.AN_flag = true;
		this.AF_flag = true;
		this.DP_flag = true;
		headerLines.stream().forEach(L->header.addMetaDataLine(L));
	}
	return header;
	}

public VariantContext apply(final VariantContext ctx) {
	if(this.recalculation_is_disabled) return ctx;
	if(!ctx.hasGenotypes()) return ctx;
	
	if(!header_was_set) throw new IllegalStateException("Vcf header was not set !");
	final Predicate<Genotype> gtFilter = G-> recalculation_ignore_gt_filtered || !G.isFiltered();
	
	final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
	if( this.DP_flag )
		{
		vcb.attribute(VCFConstants.DEPTH_KEY,
			ctx.getGenotypes().stream().
				filter(gtFilter).
				filter(G->G.hasDP()).
				mapToInt(G->G.getDP()).sum()
			);	
		}
	
	final List<Integer> acL=new ArrayList<>();
	for(final Allele alt:ctx.getAlternateAlleles())
		{
		final int ac= (int)ctx.getGenotypes().stream().
			filter(gtFilter).
			mapToLong(G->G.getAlleles().stream().filter(A->A.equals(alt)).count()).sum()
			;
		acL.add(ac);
		}
	
	final int AN= (int)ctx.getGenotypes().stream().
		filter(gtFilter).
		filter(G->G.isCalled()).
		mapToLong(G->G.getAlleles().stream().filter(A->A.isCalled()).count()).sum()
		;
	
	if(this.AC_flag) {
		vcb.rmAttribute(VCFConstants.ALLELE_COUNT_KEY);
		if(!acL.isEmpty())  vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,acL);
		}
	if(this.AN_flag) {
		vcb.rmAttribute(VCFConstants.ALLELE_NUMBER_KEY);
		vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,AN);
		}
	if(this.AF_flag) {
		vcb.rmAttribute(VCFConstants.ALLELE_FREQUENCY_KEY);
		if(AN>0 && !acL.isEmpty())
    		{
    		final List<Double> afL = new ArrayList<>(acL.size());
    		for(int x=0;x< acL.size();++x)
    			{
    			afL.add(acL.get(x)/(double)AN);
    			}
    		vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,afL);
    		}
		}	
	
	return vcb.make();
	}
}
