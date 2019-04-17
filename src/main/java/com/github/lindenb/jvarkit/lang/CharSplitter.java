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
package com.github.lindenb.jvarkit.lang;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import htsjdk.samtools.util.AbstractIterator;

/** string splitter, faster than java regex Pattern */
public interface CharSplitter {
	
	
	static class CharSplitterImpl implements CharSplitter {
		private final char delim;
		
		private CharSplitterImpl(final char delim)
			{
			this.delim = delim;
			}
		@Override
		public char getDelimiter() {
			return this.delim;
			}
		@Override
		public String toString() {
			return "CharSplitter(" + pattern() + ")";
			}
		}
/** get the delimiter */
public char getDelimiter();

/** count delimiters */
public default int countTokens(final CharSequence seq) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	int n=0,i=0;
	while(i<seqLen) {
		if(seq.charAt(i)==delim) ++n;
		i++;
		}
	return 1 + n;
	}

/** breaks the input string as a list of CharSequence */
public default List<CharSequence> splitAsCharSequenceList(final CharSequence seq) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	int prev=0;
	final List<CharSequence> list = new ArrayList<>();
	for(int i=0;i< seqLen ;i++)
		{
		if(seq.charAt(i)==delim)
			{
			list.add(seq.subSequence(prev, i));
			prev=i+1; 
			}
		}
	list.add(seq.subSequence(prev, seqLen));
	return list;
	}

/** breaks the input string as a list of CharSequence, will produce at max 'maxTokens' */
public default List<CharSequence> splitAsCharSequenceList(final CharSequence seq,final int maxTokens) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	int prev=0;
	final List<CharSequence> list = new ArrayList<>();
	for(int i=0;i< seqLen && list.size()+1 < maxTokens ;i++)
		{
		if(seq.charAt(i)==delim)
			{
			list.add(seq.subSequence(prev, i));
			prev=i+1; 
			}
		}
	list.add(seq.subSequence(prev, seqLen));
	return list;
	}

/** breaks the input string as a list of CharSequence */
public default List<String> splitAsStringList(final CharSequence seq) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	int prev=0;
	final List<String> list = new ArrayList<>();
	for(int i=0;i< seqLen ;i++)
		{
		if(seq.charAt(i)==delim)
			{
			list.add(seq.subSequence(prev, i).toString());
			prev=i+1; 
			}
		}
	list.add(seq.subSequence(prev, seqLen).toString());
	return list;
	}

/** breaks the input string as a list of CharSequence, will produce at max 'maxTokens' */
public default List<String> splitAsStringList(final CharSequence seq,final int maxTokens) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	int prev=0;
	final List<String> list = new ArrayList<>();
	for(int i=0;i< seqLen && list.size()+1 < maxTokens ;i++)
		{
		if(seq.charAt(i)==delim)
			{
			list.add(seq.subSequence(prev, i).toString());
			prev=i+1; 
			}
		}
	list.add(seq.subSequence(prev, seqLen).toString());
	return list;
	}

public default String[] split(final CharSequence seq) {
	final List<String> L = splitAsStringList(seq);
	return L.toArray(new String[L.size()]);
	}

public default String[] split(final CharSequence seq,final int maxTokens) {
	final List<String> L = splitAsStringList(seq,maxTokens);
	return L.toArray(new String[L.size()]);
	}



/** get an iterator over all the tokens in the input */
public default Iterator<CharSequence> charSequenceIterator(final CharSequence seq) {
	final char delim = getDelimiter();
	int seqLen = seq.length();
	while(seqLen-1 >=0 && seq.charAt(seqLen-1)==delim) {
		seqLen--;
		}
	final int end = seqLen;
	return new AbstractIterator<CharSequence>() {
		int curr = 0;
		int prev = 0;
		@Override
		protected CharSequence advance() {
			while(curr<end)
				{
				final char c = seq.charAt(curr);
				if(c==delim) {
					final CharSequence sub = seq.subSequence(prev, curr);
					++curr;
					prev=curr;
					return sub ;
					}
				++curr;
				}
			if(curr==end)
				{
				final CharSequence sub = seq.subSequence(prev, curr);
				++curr;
				return sub;
				}
			return null;
			}
		};
	}


/** get an iterator over all the tokens in the input */
public default Iterator<String> iterator(final CharSequence seq) {
	final Iterator<CharSequence> delegate = charSequenceIterator(seq);
	return new AbstractIterator<String>() {
		@Override
		protected String advance() {
			return delegate.hasNext()?delegate.next().toString():null;
			}
		};
	}
/** get an stream over all the tokens in the input */
public default Stream<String> stream(final CharSequence seq) {
	return StreamSupport.stream(
	          Spliterators.spliteratorUnknownSize(iterator(seq), Spliterator.ORDERED),
	          false);
	}


public default String pattern() {
	final char c = getDelimiter();
	switch(c)
		{
		case '\t': return "\\t";
		case '\n': return "\\n";
		case '\r': return "\\r";
		default: return String.valueOf(c);
		}
	}		

/** create a new CharSplitter from a string. Escaped character are interpreted */
public static CharSplitter of(final String s) {
	if(s.equals("\\t")) return TAB;
	if(s.equals("\\n")) return NEWLINE;
	if(s.equals("\\r")) return of('\r');
	if(s.isEmpty()) throw new IllegalArgumentException("cannot create a delimiter from an empty string.");
	if(s.length()!=1) throw new IllegalArgumentException("cannot create a delimiter from an empty string.");
	return of(s.charAt(0));
	}

/** create a new CharSplitter from a single character */
public static CharSplitter of(final char c) {
	return new CharSplitterImpl(c);
	}
public static final CharSplitter TAB = of('\t');
public static final CharSplitter COLON = of(':');
public static final CharSplitter COMMA = of(',');
public static final CharSplitter SPACE = of(' ');
public static final CharSplitter SEMICOLON = of(';');
public static final CharSplitter PIPE = of('|');
public static final CharSplitter DOT = of('.');
public static final CharSplitter NEWLINE = of('\n');
}
