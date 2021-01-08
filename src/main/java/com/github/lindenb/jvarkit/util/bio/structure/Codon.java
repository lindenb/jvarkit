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

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.Locatable;

public interface Codon extends TranscriptInterval {
	/** get Position of middle base of codon */
	public int getMiddle();
	/** return wether this codon is spliced */
	public default boolean isSpliced() {
		return !(getStart()+1==getMiddle() && getMiddle()+1==getEnd());
	}
	/** get Blocks (codon can be spliced */
	public default List<Locatable> getBlocks() {
		if(!isSpliced()) return Collections.singletonList(this);
		if(getStart()+1==getMiddle() && getMiddle()+1!=getEnd() ) {
			return Arrays.asList(
					new SimpleInterval(getContig(),getStart(),getMiddle()),
					new SimpleInterval(getContig(),getEnd(),getEnd())
					);
			}
		else if(getStart()+1!=getMiddle() && getMiddle()+1==getEnd() ) {
			return Arrays.asList(
					new SimpleInterval(getContig(),getStart(),getStart()),
					new SimpleInterval(getContig(),getMiddle(),getEnd())
					);
			}
		else
			{
			return Arrays.asList(
					new SimpleInterval(getContig(),getStart(),getStart()),
					new SimpleInterval(getContig(),getMiddle(),getMiddle()),
					new SimpleInterval(getContig(),getEnd(),getEnd())
					);
			}
		}
	
	
}
