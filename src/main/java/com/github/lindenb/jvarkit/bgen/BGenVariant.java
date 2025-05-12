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

import java.util.List;


import htsjdk.tribble.Feature;

public interface BGenVariant extends Feature {
	/** number of bits used for storage, or -1 if not available */
	public int getBitsPerProb();

	public default int getNAlleles() {
		return this.getAlleles().size();
		}
	public default String getAllele(int idx) {
		return this.getAlleles().get(idx);
		}
	
	
	public abstract int getPosition();
	@Override
	public default int getStart() {
		return getPosition();
		}
	@Override
	public default int getEnd() {
		return getPosition()+ getAlleles().stream().mapToInt(A->A.length()-1).max().orElse(0);
		}
	public abstract List<String> getAlleles();
	public abstract String getId();
	public abstract String getRsId();
	/** physical offset of the variant in the source bgen file; Or <=0L if not available */
	public abstract long getOffset();
	/** return true if the genotypes have been loaded, if true, getGenotypes() is available */
	public abstract boolean hasGenotypes();
	/** might throw an exception if the genotypes haven't been loaded for this variant */
	public abstract List<BGenGenotype> getGenotypes();
	/** get ith genotype */
	public default  BGenGenotype getGenotype(int i) {
		return getGenotypes().get(i);
		}
	/** get number of genotypes */
	public default int getNGenotypes() {
		return hasGenotypes()?getGenotypes().size():0;
		}
	/** get genotype by sample name */
	public BGenGenotype getGenotype(String sampleName);
	/** answer if the variant was read as phased. The information is only available if genotypes were read */
	public boolean isPhased();
}
