/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.Locatable;

public interface Transcript extends Locatable {
	/** get property map */
	public Map<String, String> getProperties();
	/** get associated gene */
	public Gene getGene();
	/** get transcript id */
	public String getId();
	/** return i-th exon scanning from 5' to 3', whatever strand */
	public Exon getExon(int index0);
	/** return i-th exon start scanning from 5' to 3', whatever strand */
	public int getExonStart(int index0);
	/** return i-th exon end scanning from 5' to 3', whatever strand */
	public int getExonEnd(int index0);
	/** return exons scanning from 5' to 3', whatever strand */
	public List<Exon> getExons();
    public char getStrand();
    public int getExonCount();
    
    
    public int getIntronCount();
	public List<Intron> getIntrons();
	/** return i-th intron scanning from 5' to 3', whatever strand */
	public Intron getIntron(int index0);

	
    public default boolean isNegativeStrand() { return getStrand()=='-';}
    public int getCdsStart();
    public int getCdsEnd();
    public int getTxStart();
    public int getTxEnd();
    /** returns true if cdsStart==cdsEnd */
	public default boolean isNonCoding()
		{
		return getCdsStart()==getCdsEnd();
		}
	
	/** get transcript length (cumulative exons sizes )*/
	public default int getTranscriptLength() {
		return getExons().stream().
				mapToInt(T->getLengthOnReference()).
				sum();
		}
	

	
}
