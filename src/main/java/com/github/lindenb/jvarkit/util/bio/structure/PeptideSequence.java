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
package com.github.lindenb.jvarkit.util.bio.structure;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;

/* a Peptide sequence */
public interface PeptideSequence<T extends CharSequence> extends CharSequence {

	/** return associated RNA */
	public T getCodingRNA();
	/** return associated genetic code */
	public GeneticCode getGeneticCode();
	
	@Override
	public default char charAt(int index) {
		char b1 = getCodingRNA().charAt(index*3 + 0);
		char b2 = getCodingRNA().charAt(index*3 + 1);
		char b3 = getCodingRNA().charAt(index*3 + 2);
		return this.getGeneticCode().translate(b1, b2, b3);
		}

	@Override
	default int length() {
		return getCodingRNA().length()/3;
		}
	
	public static <T extends CharSequence> PeptideSequence<T> of(T cDNA) {
		return of(cDNA,null);
	}
	
	public static <T extends CharSequence>  PeptideSequence<T> of(T cDNA,GeneticCode code) {
		if(code==null) code=GeneticCode.getStandard();
		return new PeptideSequenceImpl<T>(cDNA,code);
	}
	
	static class PeptideSequenceImpl<T extends CharSequence> extends AbstractCharSequence implements PeptideSequence<T> {
		final T cDNA;
		final GeneticCode geneticCode;
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
