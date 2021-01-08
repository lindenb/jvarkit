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

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

/**
 * Describes a BND structural variant
 * BND parser
 *
 */
public interface Breakend extends Locatable {
	
	
	/** Return the printed representation of this allele. */
	public String getDisplayString();
	
	/** Return  BDN delimiter :  '[' or ']'. */
	public char getDelimiter();
	/** the position of the breakend */
	public int getPosition();
	public String getLeftSequence();
	public String getRightSequence();
	
	/** return '+' is delimiter is '+' otherwise '-' */ 
	public default char getStrand() {
		return getDelimiter()=='['?'+':'-';
	}
	
	/** return true if strand is '+' */ 
	public default boolean isStrandPlus() {
		return getStrand()=='+';
	}
	
	/** return true if strand is '-' */ 
	public default boolean isStrandMinus() {
		return getStrand()=='-';
	}

	/** @return true if left/5' DNA Sequence is not empty */
	public default boolean isLeft() {
		return !StringUtils.isBlank(getLeftSequence());
	}
	
	/** @return true if right/3' DNA Sequence is not empty */
	public default boolean isRight() {
		return !StringUtils.isBlank(getRightSequence());
	}

	@Override
	public default int getStart() {
		return getPosition();
		}
	
	@Override
	public default int getEnd() {
		return getPosition();
		}
	
	/** yes, it can be a List, see vcf spec "Multiple mates" */
	public static List<Breakend> parse(final VariantContext ctx) {
		if(ctx==null || !ctx.isVariant() || !ctx.getAttributeAsString(VCFConstants.SVTYPE,"").equals("BND")) return Collections.emptyList();
		return ctx.getAlternateAlleles().
				stream().
				map(ALT->parse(ALT)).
				filter(O->O.isPresent()).
				map(O->O.get()).
				collect(Collectors.toList());
		}
	
	/** find the FIRST BND mappend on the same contig than ctx and return the interval between the ctx and the bnd.*/
	public static Optional<? extends Locatable> parseInterval(final VariantContext ctx) {
		return parse(ctx).
				stream().
				filter(B->B.contigsMatch(ctx)).
				map(B->new SimpleInterval(ctx.getContig(),
						Math.min(ctx.getStart(),B.getStart()),
						Math.max(ctx.getEnd(),B.getEnd()))).
				findFirst();
		}

	
	public static Optional<Breakend> parse(final Allele alt) {
		if(alt==null) return Optional.empty();
		if(alt.isReference()) return Optional.empty();
		if(!alt.isBreakpoint()) return Optional.empty();
		return parse(alt.getDisplayString());
	}
	
	public static Optional<Breakend> parse(final String alt) {
		if(StringUtils.isBlank(alt)) return Optional.empty();
		char delim;
		int x1 = alt.indexOf('[');
		if(x1==-1)
			{
			x1 = alt.indexOf(']');
			if(x1==-1) return Optional.empty();
			delim=']';
			}
		else
			{
			delim='[';
			}
		final int x2 = alt.indexOf(delim,x1+1);
		if(x2==-1) return Optional.empty();
		
		final int colon = alt.indexOf(':',x1+1);
		if(colon==-1 || colon >= x2) return Optional.empty();
		final String contig = alt.substring(x1+1, colon);
		if(StringUtils.isBlank(alt)) return Optional.empty();
		final int pos;
		try
			{
			pos = Integer.parseInt(alt.substring(colon+1,x2));
			if(pos<1) return Optional.empty();
			}
		catch(final NumberFormatException err) {
			return Optional.empty();
			}
		
		
		
		final BndImpl b=new BndImpl();
		b.delim = delim;
		b.displayString=alt;
		b.contig=contig;
		b.pos = pos;
		b.leftDNA = alt.substring(0,x1);
		b.rightDNA = alt.substring(x2+1);
		
		return Optional.of(b);
		}
	
	static class BndImpl implements Breakend {
		private char delim;
		private String displayString;
		private String contig;
		private int pos;
		private String leftDNA="";
		private String rightDNA="";
		
		@Override
		public char getDelimiter() {
			return this.delim;
			}
		
		@Override
		public String getContig() {
			return this.contig;
			}
		
		/** the position of the breakend */
		public int getPosition() {
			return this.pos;
		}
		

		
		@Override
		public String getDisplayString() {
			return displayString;
			}
		
		@Override
		public String getLeftSequence()  { return this.leftDNA;}
		
		@Override
		public String getRightSequence()  { return this.rightDNA;}

		
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof Breakend)) return false;
			return this.getDisplayString().equals(Breakend.class.cast(obj).getDisplayString());
			}
		
		@Override
		public int hashCode() {
			return displayString.hashCode();
			}
		
		@Override
		public String toString() {
			return getDisplayString();
			}
		}
	}
