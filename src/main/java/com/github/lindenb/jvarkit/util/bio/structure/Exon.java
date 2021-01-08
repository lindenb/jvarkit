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

import java.util.Optional;


public interface Exon extends ExonOrIntron {
	
	/** get index in getExons() */
	public int getIndex();

	/** get next intron in genomic location */
	public default Optional<Intron> getNextIntron() {
		final int idx = this.getIndex(); 
		if(idx>=this.getTranscript().getIntronCount()) return Optional.empty();
		return Optional.of(getTranscript().getIntron(idx));
		}
	/** get next intron in genomic location */
	public default Optional<Intron> getPrevIntron() {
		final int idx = this.getIndex(); 
		if(idx<=0) return Optional.empty();
		return Optional.of(getTranscript().getIntron(idx-1));
		}
	
	@Override
	public default boolean isSplicingAcceptor(final int position1)
		{
		if(!contains(position1)) return false;
		if(isPositiveStrand())
			{
			if(getIndex()== 0) return false;
			return position1==getStart();
			}
		else
			{
			if(getIndex()+1 >= getTranscript().getExonCount()) return false;
			return position1==getEnd();
			}
		}

	@Override
	public default boolean isSplicingDonor(int position1)
		{
		if(!contains(position1) || this.getLengthOnReference()<3) return false;
		if(isPositiveStrand())
			{
			if(getIndex()+1 >= getTranscript().getExonCount()) return false;
			return  (position1==getEnd()-0) ||
					(position1==getEnd()-1) ||
					(position1==getEnd()-2) ;
			}
		else
			{
			if(getIndex()== 0) return false;
			return  (position1==getStart()+0) ||
					(position1==getStart()+1) ||
					(position1==getStart()+2) ;
			}
		}
	}
