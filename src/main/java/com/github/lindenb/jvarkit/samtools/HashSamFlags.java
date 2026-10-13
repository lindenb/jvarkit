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

import htsjdk.samtools.SAMFlag;

public interface HashSamFlags {
	public int getFlags();
	public default boolean isFirstInPair() {
		return isPaired() && SAMFlag.FIRST_OF_PAIR.isSet(getFlags());
		}
	public default boolean isSecondInPair() {
		return isPaired() && SAMFlag.SECOND_OF_PAIR.isSet(getFlags());
		}
	public default boolean isPaired() {
		return SAMFlag.READ_PAIRED.isSet(getFlags());
		}
	public default boolean isProperPair() {
		return SAMFlag.PROPER_PAIR.isSet(getFlags());
		}
	public default boolean isSupplementary() {
		return SAMFlag.SUPPLEMENTARY_ALIGNMENT.isSet(getFlags());
		}
	public default boolean isDuplicate() {
		return SAMFlag.DUPLICATE_READ.isSet(getFlags());
		}
	public default boolean isUnmapped() {
		return SAMFlag.READ_UNMAPPED.isSet(getFlags());
		}
	public default boolean isMapped() {
		return !isUnmapped();
		}
	public default boolean isMateUnmapped() {
		return isPaired() && SAMFlag.MATE_UNMAPPED.isSet(getFlags());
		}
	public default boolean isMateMapped() {
		return isPaired() && !(SAMFlag.MATE_UNMAPPED.isSet(getFlags()));
		}

	public default boolean isReverseStrand() {
		return isMapped() && SAMFlag.READ_REVERSE_STRAND.isSet(getFlags());
		}
	public default boolean isMateReverseStrand() {
		return isMateMapped() && SAMFlag.MATE_REVERSE_STRAND.isSet(getFlags());
		}
	public default boolean isFailingQuality() {
		return SAMFlag.READ_FAILS_VENDOR_QUALITY_CHECK.isSet(getFlags());
		}

}
