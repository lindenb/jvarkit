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
package com.github.lindenb.jvarkit.ucsc;

import java.util.Optional;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

import htsjdk.samtools.SAMSequenceDictionary;

/**
 * A cytoband in UCSC
 *
 */
public interface Cytoband extends ExtendedLocatable,HasName {
	public String getStain();
	
	public default boolean isCentromere() {
		return getStain().contains("acen");
	}
	
	/* CSS Colors https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/lib/hCytoBand.c#L43 */
	public default String getCssColor() {
		final String stain = getStain();
		if (stain.startsWith("gneg"))
		    {
		    return "lightblue";
		    }
		else if (stain.startsWith("acen"))
		    {
		    return "orange";
		    }
		else if (stain.startsWith("gpos") && stain.length()>4)
		    {
		    int percentage;
		    try {
		    	percentage = Math.max(0, Math.min(100,Integer.parseInt(stain.substring(4))));
		    	}
		    catch(final NumberFormatException err) {
		    	percentage = 100;
		    	}
		    final int g = 40 + (int)(215.0*(percentage/100.0));
		    return "rgb("+g+","+g+","+g+")";
		    }
		else if (stain.startsWith("gvar"))
		    {
		    return "slategray";
		    }
		else 
		    {
		    return "honeydew";
		    }
		}

	/** parse line chrX	130400000	133600000	q26.2	gpos25 */
	public static Cytoband of(String[] tokens) {
		return new CytobandImpl(tokens);
		}
	
	public static Cytoband of(String line) {
		return of(CharSplitter.TAB.split(line));
		}
	
	static class CytobandImpl extends SimpleInterval implements Cytoband {
		private final String name;
		private final String stain;
		public CytobandImpl(String[] tokens) {
			super(tokens[0],1+Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]));
			this.name = tokens[3];
			this.stain = tokens[4];
			}
		@Override
		public String getName() {
			return name;
			}
		@Override
		public String getStain() {
			return this.stain;
			}
		}
	
	/** return download URL for given sam dict */
	public static Optional<String> getURLForBuild(final SAMSequenceDictionary dict) {
		final String ucsc_name;
		if(SequenceDictionaryUtils.isGRCh37(dict)) {
			ucsc_name = "hg19";
			}
		else if(SequenceDictionaryUtils.isGRCh38(dict)) {
			ucsc_name = "hg38";
			}
		else
			{
			return Optional.empty();
			}
		return Optional.of("https://hgdownload.cse.ucsc.edu/goldenpath/"+ucsc_name+"/database/cytoBand.txt.gz");
		}
}
