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
package com.github.lindenb.jvarkit.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.AbstractList;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.jexl2.JexlContext;

import com.github.lindenb.jvarkit.lang.CharSplitter;

/**
 * Header of file with fast access to the column index
 */
public class FileHeader extends AbstractList<String> {
	public static final Function<String,List<String>> DEFAULT_SPLITTER = S->CharSplitter.TAB.splitAsStringList(S);
	
	private final List<String> cols;
	private final Map<String,Integer> col2idx;
	private final Function<String,List<String>> lineSplitter;
	
	/** a Map<String,String> where we can access a column using the column header */
	public interface RowMap extends Map<String,String> {
		/** return current line as List, there may be less columns that in the header (empty columns at the end ) */
		public List<String> asList();
		/** get column content by index */
		public String at(int i);
		/** convert to an object that is suitable for JEXL expression */
		public JexlContext asJexlContext();
		}
	
	private static class JexlRow implements JexlContext {
		private final RowMapImpl row;
		JexlRow(final RowMapImpl row) {
			this.row = row;
			}
		@Override
		public Object get(final String key) {
			return row.get(key);
			}
		@Override
		public boolean has(final String col) {
			return row.owner().containsKey(col);
			}
		@Override
		public void set(final String col, Object arg1) {
			throw new java.lang.UnsupportedOperationException("cannot set "+col);
			}
		@Override
		public String toString() {
			return row.toString();
			}
		}
	
	/** class used to convert a line in a file as a Map */
	private class RowMapImpl extends AbstractMap<String,String> implements RowMap {
		private final List<String> tokens;
		protected RowMapImpl(final List<String> tokens) {
			this.tokens = Collections.unmodifiableList(tokens);
			if(tokens.size()> owner().size()) throw new IllegalArgumentException("line contains "+tokens.size()+" columns but header contains "+ owner().size()+" columns");
			}
		private FileHeader owner() {return FileHeader.this;}
		/** return this as a list of tokens in the initial list */
		@Override
		public List<String> asList() {
			return this.tokens;
			}
		@Override
		public boolean containsKey(final Object key) {
			return owner().col2idx.containsKey(key);
			}
		@Override
		public String at(final int i) {
			if(i<0 || i >= this.owner().size()) throw new IllegalArgumentException("file "+owner().size()+" columns asked 0-based index  "+ i);
			return (i<this.tokens.size()?this.tokens.get(i):"");
			}
		
		@Override
		public String getOrDefault(final Object key,final String def) {
			final int i= owner().indexOf(key);
			if(i<0) return def;
			return at(i);
			}
		
		@Override
		public String get(final Object key) {
			final int i= owner().indexOf(key);
			if(i<0) throw new IllegalArgumentException("Cannot find column \""+key+"\" in "+owner().toString());
			return at(i);
			}
		@Override
		public Set<String> keySet() {
			return owner().col2idx.keySet();//unmodifiable
			}
		@Override
		public int size() {
			return owner().size();
			}
		
		@Override
		public boolean equals(final Object o) {
			if(o==this) return true;
			if(o==null || !(o instanceof RowMapImpl)) return false;
			final RowMapImpl r = RowMapImpl.class.cast(o);
			return r.tokens.equals(this.tokens);
			}
		
		@Override
		public int hashCode() {
			return this.tokens.hashCode();
			}
		@Override
		public Set<Entry<String, String>> entrySet() {
			final Set<Entry<String, String>> set = new HashSet<>(owner().size());
			for(int i=0;i< owner().size();i++) {
				set.add(new AbstractMap.SimpleEntry<String,String>(
					owner().get(i),
					at(i)
					));
				}
			return set;
			}
		@Override
		public String toString() {
			final StringBuilder sb = new StringBuilder();
			for(int i=0;i< owner().size();i++) {
				sb.append("$").append(i+1).
					append(" ").append(owner().get(i)).
					append(" : ").append(at(i)).
					append("\n");
				}
			return sb.toString();
			}
		
		@Override
		public JexlContext asJexlContext() {
			return new JexlRow(this);
			}
		}
	
	public FileHeader(final String[] cols) {
		this(Arrays.asList(cols));
		}
	public FileHeader(final List<String> cols) {
		this(cols,DEFAULT_SPLITTER);
		}
	
	public FileHeader(final String headerLine,final CharSplitter charSplitter) {
		this(headerLine,HDR->charSplitter.splitAsStringList(HDR));
		}
	
	public FileHeader(final String headerLine,final Function<String,List<String>> lineSplitter) {
		this(Objects.requireNonNull(lineSplitter,"undefined lineSplitter").apply(Objects.requireNonNull(headerLine,"File Header is null")),lineSplitter);
		}
	
	private FileHeader(final List<String> cols,final Function<String,List<String>> lineSplitter) {
		this.cols = Collections.unmodifiableList(cols);
		this.lineSplitter = (lineSplitter==null?DEFAULT_SPLITTER:lineSplitter);
		final Map<String,Integer> hash = new HashMap<>(this.cols.size());
		for(int i=0;i< this.cols.size();i++) {
			final String col = this.cols.get(i);
			if(hash.containsKey(col)) throw new IllegalArgumentException("Duplicate column "+col+" in "+String.join("; ", cols));
			hash.put(col, i);
			}
		this.col2idx = Collections.unmodifiableMap(hash);
		}
	/** assert column exist at any position , otherwise throw exception */
	public FileHeader assertColumnExists(String colContent) {
		this.getColumnIndex(colContent);
		return this;
		}

	
	public FileHeader assertColumn(String colContent,int index) {
		if(index<0 || index>=this.size()) throw new IndexOutOfBoundsException("0<"+index+"<"+size());
		if(!get(index).equals(colContent))  throw new IllegalArgumentException("col["+index+"]="+get(index)+" but expected "+colContent);
		return this;
		}
	
	@Override
	public boolean contains(Object o) {
		if(o==null || !(o instanceof String)) return false;
		return this.containsKey(String.class.cast(o));
		}
	/** return true if column is defined in header */
	public boolean containsKey(final String colName) {
		return colName!=null && this.col2idx.containsKey(colName);
		}
	
	@Override
	public int indexOf(final Object o) {
		if(!contains(o)) return -1;
		return this.col2idx.get(String.class.cast(o));
		}
	/** get column index by column name
	 * @throws IllegalArgumentException if the column doesn't exist
	 * */
	public int getColumnIndex(final String col) {
		final Integer i = this.col2idx.get(col);
		if(i==null) {
			final StringBuilder sb = new StringBuilder("Cannot find column \"").append(col).append("\". Available are:\n");
			for(int x=0;x< this.cols.size();x++) {
				sb.append("\t$").append(x+1).append(" : \"").append(cols.get(x)).append("\"\n");
				}
			throw new IllegalArgumentException(sb.toString());
			}
		return i.intValue();
		}
	@Override
	public String get(final int index) {
		return this.cols.get(index);
		}
	@Override
	public int size() {
		return this.cols.size();
		}
	
	/** return all headers as list */
	public List<String> asList() {
		return Collections.unmodifiableList(this.cols);
		}
	
	/** convert a line in the file to a Map where keys are the columns of this header  */
	public RowMap toMap(final List<String> tokens) {
		return new RowMapImpl(tokens);
		}
	/** convert a line in the file to a Map where keys are the columns of this header  */
	public RowMap toMap(final String[] tokens) {
		return toMap(Arrays.asList(tokens));
		}
	
	/** convert a line, split it using the default splitter to a Map where keys are the columns of this header  */
	public RowMap toMap(final String row) {
		return toMap(split(row));
		}
	
	/** split a line using the internal splitter  */
	public List<String> split(final String line) {
		return this.lineSplitter.apply(line);
		}
	
	public List<RowMap> readAll(final BufferedReader br) throws IOException {
		return br.lines().map(S->toMap(S)).collect(Collectors.toList());
		}
	
	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		for(int i=0;i< this.size();i++) {
			sb.append("$").append(i+1).append(" : ").append(get(i)).append("\n");
			}
		return sb.toString();
		}
	}
