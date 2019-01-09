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
package com.github.lindenb.jvarkit.util.bio;

import java.util.HashMap;
import java.util.Map;


import com.github.lindenb.jvarkit.lang.StringUtils;

public class AminoAcids {
	
	public static interface AminoAcid
		{
		public String getName();
		public String getThreeLettersCode();
		public char getOneLetterCode();
		}
	private static class AminoAcidImpl implements AminoAcid
		{
		private final String fullName;
		private final String code3;
		private final char code1;
		AminoAcidImpl(final String fullName,final String code3,final char code1) {
			this.fullName = fullName;
			this.code3 = code3;
			this.code1 = code1;
			}
		@Override public String getName() { return this.fullName;}
		@Override public String getThreeLettersCode() { return this.code3;}
		@Override public char getOneLetterCode() { return this.code1;}
		
		@Override
		public int hashCode() {
			return this.code1*31;
			}
		
		@Override
		public boolean equals(Object obj) {
			return (obj==this);
			}
		@Override
		public String toString() {
			return this.fullName;
			}
		}
	
	public static final AminoAcid Alanine = new AminoAcidImpl("Alanine","Ala",'A');
	public static final AminoAcid Arginine = new AminoAcidImpl("Arginine","Arg",'R');
	public static final AminoAcid Asparagine = new AminoAcidImpl("Asparagine","Asn",'N');
	public static final AminoAcid Aspartate = new AminoAcidImpl("Aspartate","Asp",'D');
	public static final AminoAcid Cysteine = new AminoAcidImpl("Cysteine","Cys",'C');
	public static final AminoAcid Glutamate = new AminoAcidImpl("Glutamate","Glu",'E');
	public static final AminoAcid Glutamine = new AminoAcidImpl("Glutamine","Gln",'Q');
	public static final AminoAcid Glycine = new AminoAcidImpl("Glycine","Gly",'G');
	public static final AminoAcid Histidine = new AminoAcidImpl("Histidine","His",'H');
	public static final AminoAcid Isoleucine = new AminoAcidImpl("Isoleucine","Ile",'I');
	public static final AminoAcid Leucine = new AminoAcidImpl("Leucine","Leu",'L');
	public static final AminoAcid Lysine = new AminoAcidImpl("Lysine","Lys",'K');
	public static final AminoAcid Methionine = new AminoAcidImpl("Methionine","Met",'M');
	public static final AminoAcid Phenylalanine = new AminoAcidImpl("Phenylalanine","Phe",'F');
	public static final AminoAcid Proline = new AminoAcidImpl("Proline","Pro",'P');
	public static final AminoAcid Serine = new AminoAcidImpl("Serine","Ser",'S');
	public static final AminoAcid Threonine = new AminoAcidImpl("Threonine","Thr",'T');
	public static final AminoAcid Tryptophan = new AminoAcidImpl("Tryptophan","Trp",'W');
	public static final AminoAcid Tyrosine = new AminoAcidImpl("Tyrosine","Tyr",'Y');
	public static final AminoAcid Valine = new AminoAcidImpl("Valine","Val",'V');
	public static final AminoAcid Stop = new AminoAcidImpl("Stop","Stop",'X');

	
	@SuppressWarnings("serial")
	private static final Map<String,AminoAcid> CODE3TOAA = new HashMap<String,AminoAcid>() {{{
		put(Alanine.getThreeLettersCode(),Alanine);
		put(Arginine.getThreeLettersCode(),Arginine);
		put(Asparagine.getThreeLettersCode(),Asparagine);
		put(Aspartate.getThreeLettersCode(),Aspartate);
		put(Cysteine.getThreeLettersCode(),Cysteine);
		put(Glutamate.getThreeLettersCode(),Glutamate);
		put(Glutamine.getThreeLettersCode(),Glutamine);
		put(Glycine.getThreeLettersCode(),Glycine);
		put(Histidine.getThreeLettersCode(),Histidine);
		put(Isoleucine.getThreeLettersCode(),Isoleucine);
		put(Leucine.getThreeLettersCode(),Leucine);
		put(Lysine.getThreeLettersCode(),Lysine);
		put(Methionine.getThreeLettersCode(),Methionine);
		put(Phenylalanine.getThreeLettersCode(),Phenylalanine);
		put(Proline.getThreeLettersCode(),Proline);
		put(Serine.getThreeLettersCode(),Serine);
		put(Threonine.getThreeLettersCode(),Threonine);
		put(Tryptophan.getThreeLettersCode(),Tryptophan);
		put(Tyrosine.getThreeLettersCode(),Tyrosine);
		put(Valine.getThreeLettersCode(),Valine);
		put(Stop.getThreeLettersCode(),Stop);
		}}};
	
	public static AminoAcid getAminoAcidFromThreeLettersCode(final String s) {
		if(s==null ) return null;
		if(Stop.getThreeLettersCode().equalsIgnoreCase(s)) return Stop;
		if(s.length()!=3) return null;
		return CODE3TOAA.get(StringUtils.toTitle(s));
		}
	public static AminoAcid getAminoAcidFromOneLetterCode(final char c) {
		switch(Character.toUpperCase(c)) {
			case 'A' : return Alanine;
			case 'R' : return Arginine;
			case 'N' : return Asparagine;
			case 'D' : return Aspartate;
			case 'C' : return Cysteine;
			case 'E' : return Glutamate;
			case 'Q' : return Glutamine;
			case 'G' : return Glycine;
			case 'H' : return Histidine;
			case 'I' : return Isoleucine;
			case 'L' : return Leucine;
			case 'K' : return Lysine;
			case 'M' : return Methionine;
			case 'F' : return Phenylalanine;
			case 'P' : return Proline;
			case 'S' : return Serine;
			case 'T' : return Threonine;
			case 'W' : return Tryptophan;
			case 'Y' : return Tyrosine;
			case 'V' : return Valine;
			case 'X' : return Stop;
			default: return null;
			}
		}
}
