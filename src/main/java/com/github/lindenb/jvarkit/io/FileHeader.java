/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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

import java.util.AbstractList;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Header of file with fast access to the column index
 */
public class FileHeader extends AbstractList<String> {
	private final List<String> cols;
	private final Map<String,Integer> col2idx;
	
	/** class used to convert a line in a file as a Map */
	protected class RowMap extends AbstractMap<String,String> {
		private final List<String> tokens;
		protected RowMap(final List<String> tokens) {
			this.tokens = Collections.unmodifiableList(tokens);
			if(tokens.size()> FileHeader.this.size()) throw new IllegalArgumentException("line contains "+tokens.size()+" columns but header contains "+ FileHeader.this.size()+" columns");
			}
		/** return this as a list of tokens in the initial list */
		public List<String> asList() {
			return this.tokens;
			}
		@Override
		public boolean containsKey(final Object key) {
			return FileHeader.this.col2idx.containsKey(key);
			}
		
		@Override
		public String get(final Object key) {
			final int i= FileHeader.this.indexOf(key);
			if(i<0) return null;
			if(i>=  this.tokens.size()) return "";
			return this.tokens.get(i);
			}
		@Override
		public Set<String> keySet() {
			return FileHeader.this.col2idx.keySet();//unmodifiable
			}
		@Override
		public int size() {
			return FileHeader.this.size();
			}
		
		@Override
		public boolean equals(Object o) {
			if(o==this) return true;
			if(o==null || !(o instanceof RowMap)) return false;
			final RowMap r = RowMap.class.cast(o);
			return r.tokens.equals(this.tokens);
			}
		
		@Override
		public int hashCode() {
			return this.tokens.hashCode();
			}
		@Override
		public Set<Entry<String, String>> entrySet() {
			final Set<Entry<String, String>> set = new HashSet<>(FileHeader.this.size());
			for(int i=0;i< FileHeader.this.size();i++) {
				set.add(new AbstractMap.SimpleEntry<String,String>(
					FileHeader.this.get(i),
					(i<this.tokens.size()?this.tokens.get(i):"")
					));
				}
			return set;
			}
		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			for(int i=0;i< FileHeader.this.size();i++) {
				sb.append("$").append(i+1).
					append(" ").append(FileHeader.this.get(i)).
					append(" : ").append(this.tokens.get(i)).
					append("\n");
				}
			return sb.toString();
			}
		}
	
	public FileHeader(final String[] cols) {
		this(Arrays.asList(cols));
		}
	public FileHeader(final List<String> cols) {
		this.cols = Collections.unmodifiableList(cols);
		Map<String,Integer> hash = new HashMap<>(this.cols.size());
		for(int i=0;i< this.cols.size();i++) {
			final String col = this.cols.get(i);
			if(hash.containsKey(col)) throw new IllegalArgumentException("Duplicate column "+col);
			hash.put(col, i);
			}
		this.col2idx = Collections.unmodifiableMap(hash);
		}
	@Override
	public boolean contains(Object o) {
		if(o==null || !(o instanceof String)) return false;
		return this.col2idx.containsKey(String.class.cast(o));
		}
	@Override
	public int indexOf(Object o) {
		if(!contains(o)) return -1;
		return this.col2idx.get(String.class.cast(o));
		}
	/** get column index by column name
	 * @throws IllegalArgumentException if the column doesn't exist
	 * */
	public int getColumnIndex(final String col) {
		final Integer i = this.col2idx.get(col);
		if(i==null) throw new IllegalArgumentException("Undefined column "+col);
		return i.intValue();
		}
	@Override
	public String get(int index) {
		return this.cols.get(index);
		}
	@Override
	public int size() {
		return this.cols.size();
		}
	/** convert a line in the file to a Map where keys are the columns of this header  */
	public Map<String,String> toMap(final List<String> tokens) {
		return new RowMap(tokens);
		}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for(int i=0;i< this.size();i++) {
			sb.append("$").append(i+1).append(" : ").append(get(i)).append("\n");
			}
		return sb.toString();
		}
	}
