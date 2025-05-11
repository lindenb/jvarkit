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


public interface BGenGenotype {
	public String getSample();
	public boolean isMissing();
	public boolean isPhased();
	public double[] getProbs();
	public int getPloidy();
	public default boolean isDiploid() {
		return getPloidy()==2;
	}
	
	/** find index in getProbs of the highest prob >= treshold, or -1 if no value is>=treshold*/
	public default int findHighestProbIndex(final double treshold) {
		int best_index=-1;
		final double[] probas = getProbs();
		for(int j=0;j< probas.length;++j) {
			final double p = probas[j];
			if(p >=treshold && (best_index<0 || p> probas[best_index])) {
				best_index=j;
				}
			}
		return best_index;
		}
}
