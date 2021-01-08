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
package com.github.lindenb.jvarkit.util.bio;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;

/**
 * Describe a KOZAK sequence
 * https://en.wikipedia.org/wiki/Kozak_consensus_sequence
 * part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm
 */
public class KozakSequence
	extends AbstractCharSequence
	{
	public enum Strength {
		Strong,Moderate,Weak,nil;
		}
	
	/* kozak consensus length */
	public static final int KOZAK_LENGTH=7;
	/* position of ATG in the kozak sequence */
	public static final int KOZAK_ATG=3;

	
	private final CharSequence delegate;
	private final int offset_beg;
	
	public KozakSequence(final  CharSequence delegate,final int offset_atg) {
		this.delegate = delegate;
		if(delegate==null) throw new IllegalArgumentException("delegate is null");
		this.offset_beg = offset_atg-KOZAK_ATG;
		}
	
	@Override
	public char charAt(int idx) {
		idx += this.offset_beg;
		return idx >=0 && idx <delegate.length()?Character.toUpperCase(delegate.charAt(idx)):'N';
		}
	
	/** check there is really a ATG at KOZAK_ATG */
	public boolean hasATG() {
		if(charAt(KOZAK_ATG  )!='A') return false;
		if(charAt(KOZAK_ATG+1)!='T') return false;
		if(charAt(KOZAK_ATG+2)!='G') return false;
		return true;
		}

	/** get strength for this Kozak sequence */
	public Strength getStrength() {
		final char c0 = charAt(0);
		if(c0=='N') return Strength.nil;
		final char c6 = charAt(6);
		if(c6=='N') return Strength.nil;
		if(!hasATG()) return Strength.nil;
		final boolean ok_chart_at_0 = c0 == 'A' || c0 == 'G' ;
	
		if ( ok_chart_at_0 && c6 == 'G'){
			return  Strength.Strong;
			}
		else if (ok_chart_at_0 || c6 == 'G'){
			return  Strength.Moderate;
			}
		return Strength.Weak;
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof KozakSequence)) return false;
		return super.getString().equals(KozakSequence.class.cast(obj).getString());
	}
	
	@Override
	public int length() {
		return KOZAK_LENGTH;
		}
	}
