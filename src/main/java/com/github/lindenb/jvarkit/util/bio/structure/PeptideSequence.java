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
package com.github.lindenb.jvarkit.util.bio.structure;

import com.github.lindenb.jvarkit.bio.GeneticCode;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;

/* a Peptide sequence */
public interface PeptideSequence<T extends CharSequence> extends CharSequence {

	/** return associated RNA */
	public T getCodingRNA();
	/** return associated genetic code */
	public GeneticCode getGeneticCode();
	
	@Override
	public default char charAt(final int index) {
		final char b1 = getCodingRNA().charAt(index*3 + 0);
		final char b2 = getCodingRNA().charAt(index*3 + 1);
		final char b3 = getCodingRNA().charAt(index*3 + 2);
		return this.getGeneticCode().translate(b1, b2, b3);
		}

	@Override
	default int length() {
		return getCodingRNA().length()/3;
		}
	
	/** scan the peptide, stop when it finds a stop, otherwise return this */
	public default PeptideSequence<? extends CharSequence> getPeptideUntilStop() {
		final GeneticCode code = getGeneticCode(); 
		for(int i=0;i< this.length();i++) {
			char aa = charAt(i);
			if(code.isStop(aa)) {
				if(i+1==length()) return this;
				final CharSequence subRNA= getCodingRNA().subSequence(0, i*3+2);
				return PeptideSequence.of(subRNA, code);
				}
			}
		return this;
		}
	
	
	/** create a default Peptide with RNA Length, with default genetic code */
	public static <T extends CharSequence> PeptideSequence<T> of(T cDNA) {
		return of(cDNA,null);
	}
	/** create a default Peptide with RNA Length, with default genetic code if code is null*/
	public static <T extends CharSequence>  PeptideSequence<T> of(final T cDNA,final GeneticCode code) {
		return new PeptideSequenceImpl<T>(
				cDNA,
				code==null?GeneticCode.getStandard(): code
				);
		}
	
	static class PeptideSequenceImpl<T extends CharSequence> extends AbstractCharSequence implements PeptideSequence<T> {
		private final T cDNA;
		private final GeneticCode geneticCode;
		PeptideSequenceImpl(final T cDNA,GeneticCode geneticCode) {
			this.cDNA= cDNA;
			this.geneticCode = geneticCode;
			}
		
		@Override
		public T getCodingRNA() {
			return this.cDNA	;
			}
		@Override
		public GeneticCode getGeneticCode() {
			return geneticCode;
		}
	}
}
