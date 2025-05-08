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
package com.github.lindenb.jvarkit.bgen;

import java.util.function.ToIntFunction;

public interface BGenGenotype {
	public String getSample();
	public boolean isMissing();
	public boolean isPhased();
	public double[] getProbs();
	public int getPloidy();
	public default boolean isDiploid() {
		return getPloidy()==2;
	}
	/** return the allele indexes for this genotype 
	 * 
	 * @param fun select the index of the best GT in all the probabilities (getProbs())
	 * @return the alleles indexes, an array whith length==ploidy
	 */
	public int[] getAllelesIndexes(ToIntFunction<double[]> fun);
	/** return the allele indexes for this genotype 
	 * 
	 * @param treshold return the GT for the first probability having a proba >= treshold
	 */
	public default int[] getAllelesIndexes(double treshold) {
		return getAllelesIndexes(ARRAY->{
			final double[] p_array = getProbs();
			for(int i=0;i< p_array.length;i++) {
				if(p_array[i]>=treshold) return i;
				}
			return -1;
			});
	}
}
