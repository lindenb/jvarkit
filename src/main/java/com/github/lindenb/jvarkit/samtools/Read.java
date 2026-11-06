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
package com.github.lindenb.jvarkit.samtools;

import java.util.List;

import com.github.lindenb.jvarkit.lang.HasName;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

public interface Read extends Locatable, HasName{
	public Cigar getCigar();
	public boolean isPairedRead();
	public boolean isUnmapped();
	public boolean isDuplicate();
	public boolean isSecondary();
	public boolean isSupplementary();
	public boolean isMateUnmapped();
	public int getReferenceIndex();
	public int getMateReferenceIndex();
	public String getMateReferenceName();
	public int getMateAlignmentStart();
	public int getMateAlignmentEnd();
	public default Interval getMateInterval() {
		return new Interval(
				getMateReferenceName(),
				getMateAlignmentStart(),
				getMateAlignmentEnd(),
				isMateMinusStrand(),
				getName()
				);
		}
	
	/* inverse isUnmapped() */
	public default boolean isMapped() {
		return !isUnmapped();
		}
	
	/* inverse isUnmapped() */
	public default boolean isMateMapped() {
		return !isMateUnmapped();
		}
	
	public List<Read> getOtherCanonicalAlignments();
	
	public boolean isPlusStrand();
	public boolean isMinusStrand();
	public boolean isMatePlusStrand();
	public boolean isMateMinusStrand();

	
	public int getMappingQuality();
	public boolean isFailingVendorQualityCheck();
	/** return true if unmapped ir duplicate or secondary or supplementary */
	public default boolean isIgnorable() {
		return isUnmapped() || isDuplicate() || isSecondary() || isSupplementary() || isFailingVendorQualityCheck();
		}
	/** return true if isIgnorable() && getMappingQuality() < mapq */
	public default boolean isIgnorable(int mapq) {
		return isIgnorable() && getMappingQuality() < mapq;
		}
	/** return true if paired, read mapped, mate mapped and contig()==getMateReferenceIndex() */
	public default boolean isReadAndMateMappedOnSameChromosome() {
		return isPairedRead() && isMapped() && isMateMapped() && getReferenceIndex()==getMateReferenceIndex();
		}
}
