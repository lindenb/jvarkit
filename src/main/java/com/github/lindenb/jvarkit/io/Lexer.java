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
package com.github.lindenb.jvarkit.io;

import java.io.Closeable;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.RuntimeIOException;


/** general utility for parsing */
public class Lexer implements Closeable {
	private static final int EOF=-1;
	private final Reader reader;
	private final String inputName;
	private Lexer includeLexer = null;
	private int rowNum=1;
	private int colNum=1;
	private final List<Character> stack = new ArrayList<>();

	public Lexer(final String inputName,final Reader reader)
		{
		this.inputName = inputName;
		this.reader = reader;
		}

	public void include(final Lexer child)
		{
		if(includeLexer!=null) throw new IllegalStateException();
		this.includeLexer=child;
		}
	
	public Lexer(final Reader reader)
		{
		this("<INPUT>",reader);
		}
	
	public String getLocation()
		{
		if(this.includeLexer!=null) return this.includeLexer.getLocation();
		return this.inputName+" row: "+this.rowNum+" column:"+this.colNum;
		}
	
	public String toString()
		{
		return getLocation();
		}
	
	private int _read()
		{
		try
			{
			final int c = reader.read();
			switch(c)
				{
				case EOF: break;
				case '\n': rowNum++;colNum=1;break;
				default: colNum++;break;
				}
			return c;
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}		
		}
	public int peek(final int pos)
		{
		while(pos <= this.stack.size())
			{
			final int c = this._read();
			if(c==-1) return EOF;
			switch(c)
				{
				case '\n': rowNum++;colNum=1;break;
				default: colNum++;break;
				}
	public int get(final int pos)
		{
		while(pos<=stack.size())
			{
			final int c = this._read();
			if(c==EOF) return EOF;
			this.stack.add((char)c);
			}
		return this.stack.get(pos);
		}

	public int peek() { return peek(0); }
	public int consumme(int n) {
		int c=0;
		while(n>0 && !stack.isEmpty())
			{
			this.stack.remove(0);
			n--;
			c++;
			}
		while(n>0 && this._read()!=EOF)
			{
			n--;
			c++;
			}
		return c;
		}
	
	public boolean eof() { return get()!= EOF;}
	
	public int skip(final Predicate<Character> filter) {
		int c,n=0;
		while((c=get())!=EOF && filter.test((char)c))
			{
			consumme();
			++n;
			}
		return n;
		}
	public int skipws() {return skip(C->Character.isWhitespace(C));}
	
	public String peekLine()
		{
		StringBuilder sb=null;
		int c;
		int i=0;
		while((c=get(i))!= EOF && c!='\n')
			{
			if(sb==null) sb=new StringBuilder();
			sb.append((char)c);
			i++;
			}
		return sb==null?null:sb.toString();
		}
	
	public String readLine()
		{
		StringBuilder sb=null;
		int c;
		while((c=get())!= EOF && c!='\n')
			{
			if(sb==null) sb=new StringBuilder();
			sb.append((char)c);
			consumme(1);
			}
		return sb==null?null:sb.toString();
		}
	
	public String matcher(final Pattern p)
		{
		final StringBuilder sb= new StringBuilder();
		Matcher matcher= p.matcher("");
		int c;
		while((c=get())!=-1 && c!='\n')
			{
			sb.append((char)c);
			matcher.reset(sb.toString());
			if(!matcher.matches()) break;
			}
		
		if(matcher.hitEnd())
			{
			return sb.toString();
			}	
		return null;
		}

	public boolean downstream(final int pos,final String s) {
	for(int i=0;i< s.length();i++) {
		final int c = peek(pos+i);
		if(c==EOF ||c!=(int)s.charAt(i)) return false;
	}
	return true;
}

	public boolean downstream(final String s) {
		Objects.requireNonNull(s);
		return downstream(0, s);
		}

	
	@Override
	public void close() throws IOException {
		this.reader.close();
		}
	}
