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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.samtools.util.CoordMath;

/**
 * A class used to compare SV on their size
 * The fraction of overlap is not the same when comparing the large and the small SV.
 */
public class OverlapComparator  {
	public final static String OPT_DESC="Two SV have are the same if they share a fraction 'x' of their bases. For very small SV the fraction can be quite small while for large SV the fraction should be close to 1. "+
			"The Syntax is the following : (<MAX_SIZE_INCLUSIVE>:<FRACTION as double or percent>;)+ . For example if the SV as a size of 99bp, the fraction used with be 0.6 for '10:0.1;100:0.6;1000:0.9'. " +
			"For the smallest size, a simple overlap is a positive match.";
	public final static String DEFAULT_VALUE = "10:0.5;100:75%;1000:80%;10000:90%";

	private static final class Pair {
		int size;
		double fraction;
		}
	public final List<Pair> conditions = new ArrayList<>();
	
	/** create a default OverlapComparator */
	public static OverlapComparator makeDefault() {
		return parse(DEFAULT_VALUE);
		}
	
	/** parse a OverlapComparator from String */
	public static OverlapComparator parse(final String query) {
		try {
			final OverlapComparator oc = new OverlapComparator();
			final FractionConverter fc= new FractionConverter();
			for(String s: CharSplitter.SEMICOLON.split(query)) {
				final String[] s2 = CharSplitter.COLON.split(s,2);
				if(s2.length!=2) throw new IllegalArgumentException("colon missing in "+query);
				final Pair p=new Pair();
				p.size = Integer.parseInt(s2[0]);
				if(p.size<=0) throw new IllegalArgumentException("length should be > 0 in "+query);
				if(oc.conditions.stream().anyMatch(C->C.size==p.size)) throw new IllegalArgumentException("duplicate size "+p.size+" in: "+query);
				p.fraction = fc.applyAsDouble(s2[1]);
				oc.conditions.add(p);
				}
			if(oc.conditions.isEmpty()) throw new IllegalArgumentException("No definition was found in :"+query);
			Collections.sort(oc.conditions,(A,B)->Integer.compare(A.size, B.size));
			return oc;
			}
		catch(Throwable err) {
			throw new IllegalArgumentException("Cannot parse OverlapComparator : "+query,err);
			}
		}
	
	/** String Converter for JCommander */
	public static class StringConverter implements IStringConverter<OverlapComparator> {
		@Override
		public OverlapComparator convert(final String query) {
			return OverlapComparator.parse(query);
			}
		}
	
	private OverlapComparator() {
		}
	
	private Pair findPairBySize(int len) {
		for(int i= this.conditions.size()-1;i>=0;i--) {
			final Pair p= this.conditions.get(i);
			if(len>=p.size) return p;
			}
		return null;
		}
	
	// I don't use locatable because the contig might not match chr1 vs '1' 
	/** test of two coordinates overlap */
	public boolean test(int start1,int end1, int start2, int end2) {
		if(!CoordMath.overlaps(start1, end1, start2, end2)) return false;
		final int len1= CoordMath.getLength(start1, end1);
		final int len2= CoordMath.getLength(start2, end2);
		final int beg  = Math.max(start1, start2);
		final int stop = Math.min(end1, end2);
		final double comm = CoordMath.getLength(beg,stop);
		
		
		final double f1 = comm/len1;
		Pair p = findPairBySize(len1);
		if(p!=null && f1< p.fraction) return false;
		
		final double f2 = comm/len2;
		p = findPairBySize(len2);
		if(p!=null && f2 < p.fraction) return false;
		return true;
		}
	
	@Override
	public String toString() {
		return conditions.stream().map(C->String.valueOf(C.size)+":"+C.fraction).collect(Collectors.joining(";"));
		}

}
