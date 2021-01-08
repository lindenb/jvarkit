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
package com.github.lindenb.jvarkit.hic;

import java.io.Closeable;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;


/*
<pre>
## Header

magic	String	'HIC'
version	int
master-position long --> file-offset to FOOTER
genomeId  String
number-of-attributes  int
   attribute-key String
   attribute-value String
number-of-contigs  int
   contig-name String
   contig-length int
nBP int number of base pair resolutions
   resolution int
NF int number of fragment resolutions 
   resolution int


## Footer


## BLOCK at (position,size) version '8'

   read the array at 'position' and inflate (it is compressed)

n-items in number of items in this block
binXOffset int
binYOffset int
useShort byte  
type byte (1 or 2)






</pre>

*/

/** interface for a Hi-C reader */
public interface HicReader extends Closeable {
	
/** get source of this reader (path, url...) or null */
public Object getSource();	

/** get dictionary */
public SAMSequenceDictionary getDictionary();

/** get genome build */
public String getBuild();

/** get attributes */
public Map<String,String> getAttributes();

/** get version of hic format */
public int getVersion();

/** get the base pair resolutions */
public Set<Integer> getBasePairResolutions();

/** get the fragment resolutions */
public Set<Integer> getFragmentResolutions();

/** try to parse an interval as string */
public Optional<Locatable> parseInterval(final String s);

public interface QueryCallBack {
	
	public void reportContact(
			String contig1,int start1,int end1,
			String contig2,int start2,int end2,
			final Normalization norm,
			final Unit unit,
			final int binsize, 
			final float value
			);
 	
	default void warning(Object o) {
		System.err.println("[WARN]"+o);
		}
	default void error(Object o) {
		System.err.println("[ERROR]"+o);
		}
}

/** query the hic file */
public boolean query(
		final Locatable interval1,
		final Locatable interval2,
		final Normalization norm,
		final int binsize, 
		final Unit unit,
		final QueryCallBack callback
		);
}



