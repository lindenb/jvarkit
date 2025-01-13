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
package com.github.lindenb.jvarkit.util.bio.gtf;

import htsjdk.tribble.Feature;

import java.util.Map;

import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;

public interface GTFLine 
	extends ExtendedLocatable,Feature

	{
	public static final int NO_PHASE=-1;
	public static final char NO_STRAND='.';
	
	/** get the original line of the GTF/GFF line */
	public String getLine();
	public String getSource();
	public String getType();
	public Double getScore();
	public char getStrand();
	public Map<String, String> getAttributes();
	public String getAttribute(final String key);
	public int getPhase();
	
	public default boolean hasAttribute(final String key) {
		return getAttributes().containsKey(key);
		}
	public default String getGeneName() {
		return getAttribute("gene_name");
		}
	public default String getGeneId() {
		return getAttribute("gene_id");
		}
	public default String getTranscriptId() {
		return getAttribute("transcript_id");
		}
	public default boolean isGene() {
		return getType().equals("gene");
		}
	public default boolean isTranscript() {
		return getType().equals("transcript");
		}
	public default boolean isExon() {
		return getType().equals("exon");
		}
	public default boolean isCDS() {
		return getType().equals("CDS");
		}
	public default boolean isPostiveStrand() {
		return getStrand()=='+';
		}
	public default boolean isNegativeStrand() {
		return getStrand()=='-';
		}
	}
