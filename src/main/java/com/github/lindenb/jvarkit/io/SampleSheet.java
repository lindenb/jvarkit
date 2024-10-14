/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.RuntimeIOException;


public interface SampleSheet extends List<FileHeader.RowMap>{
	/** get file header */
	public FileHeader getHeader();
	/** get source of samplesheet */
	public String getSource();
	/** return all values for a given table */
	public default List<String> getColumnContentByName(final String colName) {
		getHeader().assertColumnExists(colName);
		final List<String> L = new ArrayList<>(this.size());
		for(int i=0;i< this.size();i++) {
			L.add(this.get(i).get(colName));
			}
		return L;
		}
	/** return all values for a given table */
	public default Set<String> getColumnContentByNameAsSet(final String colName) {
		getHeader().assertColumnExists(colName);
		final Set<String> set = new HashSet<>();
		for(int i=0;i< this.size();i++) {
			set.add(this.get(i).get(colName));
			}
		return set;
		}
	/** assert column has unique values . Return this.*/
	public default SampleSheet assertColumnHasUniqueValues(final String colName) {
		getHeader().assertColumnExists(colName);
		final Map<String,Integer> s2row = new HashMap<>(this.size());
		for(int i=0;i< size();i++) {
			final String s=this.get(i).get(colName);
			if(s2row.containsKey(s)) throw new IllegalStateException("assertion failed. Value \""+s+"\" was found twice in row["+i+"] and row["+s2row.get(s)+"] in "+getSource());
			s2row.put(s, i);
			}
		return this;
		}
	public default void writeMarkdown(Appendable out) {
		try {
		final int[] col2size = new int[this.getHeader().size()];
			for(int x=0;x< getHeader().size();x++) {
				col2size[x] = getHeader().get(x).length();
				}
			for(int i=0;i< size();i++) {
				for(int x=0;x< getHeader().size();x++) {
					col2size[x] = Math.max(col2size[i],this.get(i).at(x).length());
					}
				}
			out.append("| ");
			for(int x=0;x< getHeader().size();x++) {
				if(x>0) out.append(" | ");
				String format = "%-"+col2size[x]+"s";
				out.append(String.format(format, getHeader().get(x)));
				}
			out.append(" |\n");
			
			out.append("|");
			for(int x=0;x< getHeader().size();x++) {
				if(x>0) out.append("|");
				out.append(StringUtils.repeat(col2size[x]+2, '-'));
				}
			out.append("|\n");
			
			
			out.append("| ");
			for(int i=0;i< size();i++) {
				
				for(int x=0;x< getHeader().size();x++) {
					if(x>0) out.append(" | ");
					String format = "%-"+col2size[x]+"s";
					out.append(String.format(format, get(i).at(x)));
					}
				}
			out.append("|\n");
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	}
