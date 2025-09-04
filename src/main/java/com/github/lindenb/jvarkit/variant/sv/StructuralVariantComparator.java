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
package com.github.lindenb.jvarkit.variant.sv;

import java.util.List;
import java.util.function.BiPredicate;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.locatable.SimplePosition;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

/** test if two SV are the same */
public class StructuralVariantComparator /* implements BiPredicate<VariantContext,VariantContext> see https://github.com/cbeust/jcommander/issues/582 */{

	@Parameter(names={"--bnd-distance"},description="Two BND variants are the same if their bounds are distant by less than xxx bases. "+ DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class ,splitter=com.github.lindenb.jvarkit.jcommander.NoSplitter.class)
	private int bnd_max_distance = 100;
	@Parameter(names={"--sv-fraction"},description= OverlapComparator.OPT_DESC,splitter = NoSplitter.class,converter = OverlapComparator.StringConverter.class)
	private OverlapComparator overlap_comparator = OverlapComparator.makeDefault();
	@Parameter(names={"--force-svtype"},description="When comparing two SV variants, their INFO/SVTYPE should be the same. Default is to just use coordinates to compare non-BND variants.")
	private boolean svtype_mus_be_the_same = false;
	@Parameter(names={"--sv-alleles-bases"},description="When comparing two non-BND SV variants, use their ALT alleles to adjust the interval. It solves the problem of  'chr2:10556028:AATTATATAATTAAATTAATTATATAATT:A'  vs 'chr2:10556028:A:AATTATATAATTAAATTAATTATATAATT'. See https://twitter.com/yokofakun/status/1169182491606999046 ")
	private boolean use_bases = false;
	@Parameter(names={"--check-bnd-mate"},description="When comparing two BND, check that their mate (using the ALT allele) are the same too")
	private boolean check_bnd_mate = false;

	/** test if two contigs are the same */
	private BiPredicate<String,String> contigComparator = (A,B)->A.equals(B);
	
public StructuralVariantComparator() {
		
	}

public StructuralVariantComparator setContigComparator(BiPredicate<String,String> contigComparator) {
	this.contigComparator = contigComparator;
	return this;
	}

/** must be accessible when using vcf.query(... ) */
public int getBndDistance() {
	return bnd_max_distance;
	}

private String getSvtype(final VariantContext ctx) {
	return ctx.getAttributeAsString(VCFConstants.SVTYPE, VCFConstants.MISSING_VALUE_v4);
	}
/** length of allele is it's a DNA or -1 */
private int getATGCNLength(final Allele a) {
	if(a.isSymbolic() || a.isNoCall() || a.equals(Allele.SPAN_DEL) ) return -1;
	final String s = a.getDisplayString();
	if(StringUtils.isBlank(s)) return -1;
	if(!AcidNucleics.isATGCN(s)) return -1;
	return s.length();
	}


/** convert variant to interval, using DNA bases if needed  */
private SimpleInterval toInterval(final VariantContext ctx) {

	
	if(this.use_bases && ctx.getNAlleles()>1) {
		final int x1 = ctx.getStart();
		int x2 = ctx.getEnd();
		for(final Allele a: ctx.getAlleles()) {
			final int len = getATGCNLength(a);
			if(len == -1) return  new  SimpleInterval(ctx);
			final int x3 = x1 + len -1;
			x2 = Math.max(x2, x3);
			}
		return new SimpleInterval(ctx.getContig(),x1,x2);
		}
	else
		{
		return new  SimpleInterval(ctx);
		}
	}

private boolean withinDistanceOf(final Locatable a,final Locatable b, final int distance) {
    return this.contigsMatch(a,b) && 
            CoordMath.overlaps(a.getStart(), a.getEnd(), b.getStart()-distance, b.getEnd()+distance);
	}




private boolean testBNDDistance(final Locatable a,final Locatable b) {
	final Function<Locatable,Locatable> bnd2interval=LOC->{
		if(!(LOC instanceof VariantContext)) {
			return new SimplePosition(LOC.getContig(), LOC.getStart());
			}
		final VariantContext CTX = VariantContext.class.cast(LOC);
		int extend_5 = 0;//can be negative
		int extend_3 = 0;//can be negative
		if(CTX.hasAttribute("CIPOS")) {
			try {
				// can be negative
				final List<Integer> xList=CTX.getAttributeAsIntList("CIPOS", 0);
				if(xList.size()==2) {
					extend_5 = xList.get(0).intValue();
					extend_3 = xList.get(1).intValue();
					if(extend_5>0) extend_5=0;
					if(extend_3<0) extend_3=0;
					}
				}
			catch(final Throwable err)
				{
				extend_5=0;
				extend_3=0;
				}
			}
		return new SimpleInterval(
				CTX.getContig(),
				CTX.getStart()+extend_5/* can be negative */,
				CTX.getEnd()+extend_3
				);
		};
	
	//  a.withinDistanceOf(b, this.small_length_on_ref ); <-- NO contigs can be renamed
	// we use SimplePosition because getEnd() can be large when the BND is on the same chromosome
	return this.withinDistanceOf(
			bnd2interval.apply(a),
			bnd2interval.apply(b), 
			this.bnd_max_distance);
	}

/** test if two locatable have overlaping overlaps after contig name conversion. This is NOT the function used to compare two SV.*/
public boolean testSimpleOverlap(final Locatable a,final Locatable b) {
	return this.withinDistanceOf(a, b,0);
	}


/** return true if SV types will be tested to be the same */
public boolean isTestingSvTypes() {
	return this.svtype_mus_be_the_same;
	}

/** return true two contigs are the same, after conversion */
private boolean contigsMatch(final Locatable a,final Locatable b)  {
	if(a.contigsMatch(b)) return true;
	return this.contigComparator.test(a.getContig(), b.getContig());
}

// https://github.com/cbeust/jcommander/issues/582
public BiPredicate<VariantContext,VariantContext> asBiPredicate() {
	return (V1,V2)->this.test(V1, V2);
}


//@Override
public boolean test(final VariantContext a, final VariantContext b) {
	if(a==null || b==null) return false;
	
	if(!this.contigsMatch(a,b)) return false;

	final String typea = getSvtype(a);
	if(typea.equals(VCFConstants.MISSING_VALUE_v4)) return false;
	final String typeb = getSvtype(b);
	if(typeb.equals(VCFConstants.MISSING_VALUE_v4)) return false;

	
	if(isTestingSvTypes() && !typea.equals(typeb)) return false;
	

	if(typea.equals("BND") &&  typeb.equals("BND")) {
		if(! testBNDDistance(a,b)) {
			return false;
			}
		if(this.check_bnd_mate) {
			final List<Breakend> bk1 =  Breakend.parse(a);
			final List<Breakend> bk2 =  Breakend.parse(b);
			if(!bk1.isEmpty() && ! bk2.isEmpty())
				{
				return bk1.stream().anyMatch(BK1->bk2.stream().anyMatch(BK2->testBNDDistance(BK1, BK2)));
				}
			
			final String chr2a = a.getAttributeAsString("CHR2", "");
			final String chr2b = b.getAttributeAsString("CHR2", "");
			if(!StringUtils.isBlank(chr2a) && !StringUtils.isBlank(chr2b)) {
				if(!((chr2a.equals(chr2b) || this.contigComparator.test(chr2a, chr2b)))) {
					return false;
					}
				//TODO
				}
			return false;
			}
		
		return true;
		}
	else if(typea.equals("BND") ||  typeb.equals("BND"))
		{
		return false;
		}
	else
		{
		final SimpleInterval interval1 =toInterval(a);
		final SimpleInterval interval2 = toInterval(b);
		if(!testSimpleOverlap(interval1,interval2)) return false;

		
		return overlap_comparator.test(
				interval1.getStart(),interval1.getEnd(),
				interval2.getStart(),interval2.getEnd()
				);
		}
	}
}
