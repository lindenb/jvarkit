/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.par;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

/** find human pseudo autosomal regions PAR */
public interface PseudoAutosomalRegion {
public enum Label {autosomal,pseudoautosomal,sexual,mixed /* overlap PAR and sex */,mitochondrial,undefined}
public static Optional<PseudoAutosomalRegion> getInstance(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return Optional.empty();
	ParRegionDefition parDefinition = null;
	// https://en.wikipedia.org/wiki/Pseudoautosomal_region#Location
	if(SequenceDictionaryUtils.isGRCh38(dict)) {
		parDefinition = new ParRegionDefition(
			"GRCh38",
			new int[] {10_001,2_781_479, 155_701_383,156_030_895},
			new int[] {10_001,2_781_479, 56_887_903, 57_217_415}
			);
		}
	else if(SequenceDictionaryUtils.isGRCh37(dict)) {
		parDefinition = new ParRegionDefition(
			"GRCh37",
			new int[] {60_001,2_699_520, 154_931_044,155_260_560 },
			new int[] {10_001,2_649_520, 59_034_050,59_363_566}
			);
		}
	return Optional.ofNullable(parDefinition);
	}

public Label getLabel(final Locatable l);
public String getDescription();
public String getName();
/** split src between region that are 100% autosomal or 100% sexual */ 
public List<Locatable> split(final Locatable src);
 
static class ParRegionDefition implements PseudoAutosomalRegion
	{
	private final String name;
	private final int[] Xrgn;
	private final int[] Yrgn;
	private final Pattern autosomalMatcher= Pattern.compile("(chr)?[0-9]+",Pattern.CASE_INSENSITIVE);
	private final Pattern mitochondrialMatcher= Pattern.compile("(chr)?m(t)?+",Pattern.CASE_INSENSITIVE);
	
	ParRegionDefition(final String name,int[] Xrgn, int[] Yrgn) {
		if(Xrgn.length%2!=0 || Xrgn.length != Yrgn.length) throw new IllegalArgumentException();
		this.name = name;
		this.Xrgn = Xrgn;
		this.Yrgn = Yrgn;
		}
	
	@Override
	public int hashCode() {
		return name.hashCode();
		}
	
	public Label getLabel(final Locatable src) {
		if(src==null) throw new IllegalArgumentException("src is null");
		final int[] array;
		if(isX(src.getContig())) {
			array = this.Xrgn;
			}
		else if(isY(src.getContig())) {
			array = this.Yrgn;
			}
		else if(src.getContig().startsWith("chrX_") || src.getContig().startsWith("chrY_")) {
			return Label.sexual;
			}
		else if(mitochondrialMatcher.matcher(src.getContig()).matches()) {
			return Label.mitochondrial;
			}
		else if(!autosomalMatcher.matcher(src.getContig()).matches()) {
			return Label.undefined;
			}
		else
			{
			return Label.autosomal;
			}
		for(int i=0;i+1 <array.length;i+=2) {
			final int parStart  = array[i];
			final int parEnd  = array[i+1];
				if(CoordMath.overlaps(parStart, parEnd, src.getStart(), src.getEnd())) {
					if(CoordMath.encloses(parStart, parEnd, src.getStart(), src.getEnd())) {
						return Label.pseudoautosomal;
						}
					return Label.mixed;
				}	
			}
		return Label.sexual;
		}
	private boolean isX(final String contig) {
		return contig.equals("chrX") || contig.equals("X");
		}
	
	private boolean isY(final String contig) {
		return contig.equals("chrY") || contig.equals("Y");
		}
	private List<Locatable> split2(final Locatable src,final int parStart,int parEnd) {
		if(!CoordMath.overlaps(parStart, parEnd, src.getStart(), src.getEnd()) ||
			CoordMath.encloses(parStart, parEnd, src.getStart(), src.getEnd())) {
			return Collections.singletonList(src);
			}
		final int x2 = Math.max(parStart, src.getStart());
		final int x3 = Math.min( src.getEnd(),parEnd);
		final List<Locatable> L = new ArrayList<>();
		Locatable l;
		if(src.getStart() < x2-1) {
			l = new SimpleInterval(src.getContig(),src.getStart(),x2-1);
			L.add(l);
			}
		l = new SimpleInterval(src.getContig(),x2,x3);
		L.add(l);
		if(x3+1 < src.getEnd()) {
			l = new SimpleInterval(src.getContig(),x3+1,src.getEnd());
			}
		L.add(l);
		
		return L;
		}
	@Override
	public List<Locatable> split(final Locatable src) {
		final int[] array;
		if(isX(src.getContig())) {
			array = this.Xrgn;
			}
		else if(isY(src.getContig())) {
			array = this.Yrgn;
			}
		else
			{
			return Collections.singletonList(src);
			}
		List<Locatable> l = Collections.singletonList(src);
		for(int i=0;i+1 <array.length;i+=2) {
			final List<Locatable> l2 = new ArrayList<>();
			for(Locatable src2:l) {
				l2.addAll(split2(src2,array[i],array[i+1]));
				}
			l=l2;
			}
		return l;
		}
	
	private String intervals(final String ctg,int[] array) {
		String s="";
		for(int i=0;i+1 <array.length;i+=2) {
			s+=ctg+":"+array[i]+"-"+array[i+1]+" ";
			}
		return s;
		}
	@Override
	public String getDescription() {
		return "Interval is mapped on a sexual chromosome on "+getName()+" excluding PAR region(s) on "+intervals("X",Xrgn)+intervals("Y",Yrgn); 
		}
	@Override
	public String getName() {
		return this.name;
		}
	@Override
	public String toString() {
		return getName();
		}
	}
}
