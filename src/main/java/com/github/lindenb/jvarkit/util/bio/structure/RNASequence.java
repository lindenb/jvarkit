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

import java.util.OptionalInt;


/** base interface for the RNA of a transcript */
public interface RNASequence extends CharSequence,StrandedLocatable {
/** get Associated Transcript */
public Transcript getTranscript();
@Override
public default String getContig() {
		return getTranscript().getContig();
	}
/** returns start of getTranscript*/
@Override
public default int getStart() {
		return getTranscript().getStart();
	}

/** returns end of getTranscript*/
@Override
public default int getEnd() {
		return getTranscript().getEnd();
	}
@Override
public default char getStrand() {
	return getTranscript().getStrand();
}
/** convert 0-based position in genomic to position 0-based in RNA. */
public OptionalInt convertGenomic0ToRnaIndex0(int g0);

/** convert 0-based position in RNA to position 0-based in genomic. */
public int convertRnaIndex0ToGenomic0(int rna0);

}
