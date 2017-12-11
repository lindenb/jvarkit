/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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


import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;

/**
 * GFF/GTF Codec
 * @author lindenb
 *
 */
public class GTFCodec extends AsciiFeatureCodec<GTFLine>{
	private static final Pattern tab=Pattern.compile("[\t]");
	private static final String GFF_VERSION="##gff-version";
	private GTFHeaderImpl header=null;

	public static interface GTFHeader
		{
		public boolean isGff3();
		public List<String> getLines();
		}

	
	public static class GTFHeaderImpl implements GTFHeader
		{
		private final List<String> lines = new ArrayList<>();
		@Override
		public boolean isGff3() {
				for(final String line:this.lines) {
				if(line.startsWith(GFF_VERSION+" "))
					{
					final String version =line.substring(GFF_VERSION.length()).trim();
					if(version.equals("3"))
						{
						return true;
						}
					}
				}
			return false;
			}
		@Override
		public List<String> getLines() {
			return this.lines;
			}
		@Override
		public String toString() {
			return String.join("\n", this.lines);
			}
		}

	
	public GTFCodec() {
		super(GTFLine.class);
		}
	
	@Override
	public boolean canDecode(final String path) {
		if(StringUtil.isBlank(path)) return false;
		return true;
		}

	@Override
	public GTFLine decode(final LineIterator r)
		{
		for(;;)
			{
			if(!r.hasNext()) return null;
			final String line=r.next();
			if(line.startsWith("#")) continue;
			final GTFLine record =  decode(line);	
			if(record==null) continue;
			return record;
			}
		}
		
	@Override
	public GTFHeader readActualHeader(final LineIterator r) {
		if(this.header!=null) throw new RuntimeIOException("Reader already read");
		this.header = new GTFHeaderImpl();
		while(r.hasNext() && r.peek().startsWith("#"))
			{
			
			this.header.lines.add(r.next());
			}
		return this.header;
		}
		
		
		
	public GTFLine decode(final String line)
		{
		/* non, on s'en fout à vrai dire...
		if(this.header==null) {
			throw new RuntimeIOException("header was not parsed");
		}*/
		if(line.startsWith("#") || line.isEmpty()) return null;
		return new GTFLineImpl(GTFCodec.tab.split(line));
		}
	
	
	private static class GTFLineImpl implements GTFLine
		{
		final String tokens[];
		final int start;
		final int end;

		public GTFLineImpl(final String tokens[])
			{
			this.tokens = tokens;
			if(tokens.length<8)
				{	
				throw new JvarkitException.TokenErrors("Expected 8 columns",tokens);
				}
			this.start = Integer.parseInt(tokens[3]);
			this.end = Integer.parseInt(tokens[4]);
			}
		
		@Override
		public int hashCode() {
			return Arrays.hashCode(this.tokens);
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof GTFLine)) return false;
			return this.getLine().equals(GTFLine.class.cast(obj).getLine());
		}
		
		private String get(int col) {
			return this.tokens[col];
		}
		
		public String getContig() {
			return get(0);
		}
		
		public String getSource() {
			return get(1);
		}
		
		public String getType() {
			return get(2);
		}
		
		public int getStart() {
			return start;
		}
		
		public int getEnd() {
			return end;
		}
	
		public Double getScore() {
			return get(5).equals(".")?null:Double.parseDouble(get(5));
		}
	
		public char getStrand() {
			return get(6).charAt(0);
		}
	
	
		public int getPhase() {
			return (get(7).equals(".")?GTFLine.NO_PHASE:Integer.parseInt(get(7)));
		}
		
		@Override
		public String getLine() {
			return String.join("\t", this.tokens);
			}
		@Override
		public Iterator<Entry<String, String>> iterator() {
			return new AttIter(get(8));
			}
		
		@Override
		public String getAttribute(final String key) {
			for(final Iterator<Map.Entry<String,String>> iter=this.iterator();iter.hasNext();)
				{
				final Map.Entry<String,String> kv = iter.next();
				if(kv.getKey().equals(key)) return kv.getValue();
				}
			return null;
			}
		
		@Override
		public Map<String, String> getAttributes() {
			final Map<String,String> hash = new LinkedHashMap<>();
			for(final Iterator<Map.Entry<String,String>> iter=this.iterator();iter.hasNext();)
				{
				final Map.Entry<String,String> kv = iter.next();
				hash.put(kv.getKey(), kv.getValue());
				}
			return hash;
			}
		
		@Override
		public String toString() {
			return getLine();
			}
		}

	private static class AttIter extends AbstractIterator<Map.Entry<String, String>> {
		private final String mapStr;
		private int k=0;
		AttIter(final String mapStr)
			{
			this.mapStr = mapStr;
			}
		private void skipws() {
			while( this.k < this.mapStr.length() &&
				Character.isWhitespace(this.mapStr.charAt(this.k)))
				{
				++this.k;
				}
			}
		
		@Override
		protected Entry<String, String> advance() {
			skipws();
			if(k>=this.mapStr.length()) return null;
			for(;;)
				{
				skipws();
				if(this.k>=this.mapStr.length()) return null;
				char c= mapStr.charAt(k);
				if(c==';') { ++k; continue;}
				/* read KEY */
				final StringBuilder sbk=new StringBuilder();
				while( this.k < mapStr.length()) {
					c= mapStr.charAt(k);
					++k;
					if(c=='=' || Character.isWhitespace(c))
						{
						break;
						}
					sbk.append(c);
					}
				/* SKIP WS */
				skipws();
				/* EQUAL SIGN */
				if( this.k < mapStr.length() && mapStr.charAt(k)=='=') {
					++k;
					}
				/* SKIP WS */
				skipws();
				
				if( this.k >= mapStr.length())
					{
					if(sbk.length()==0) return null;
					return new AbstractMap.SimpleEntry<String,String>(sbk.toString(),"");
					}
				
				/* read VALUE */
				final StringBuilder sbv=new StringBuilder();
				c=(this.k < mapStr.length()?this.mapStr.charAt(this.k):'\0');
				// quoted string
				if( c == '\"')
					{
					++this.k;
					while( this.k < mapStr.length()) {
						c= mapStr.charAt(k);
						++k;
						if(c=='\\')
							{
							c=(k < mapStr.length()?mapStr.charAt(k):'\0');
							++k;
							switch(c) {
								case '"': sbv.append("\"");break;
								case '\'': sbv.append("\'");break;
								case 't': sbv.append("\t");break;
								case 'n': sbv.append("\n");break;
								default:break;
								}
							}
						else if(c=='\"')
							{
							break;
							}
						else
							{
							sbv.append(c);
							}
						}
					}
				else
					{
					while( this.k < this.mapStr.length()) {
						c= this.mapStr.charAt(k);
						++k;
						if(c==';' || Character.isWhitespace(c))
							{
							break;
							}
						sbv.append(c);
						}
					}
				final AbstractMap.SimpleEntry<String,String> entry= new AbstractMap.SimpleEntry<String,String>(sbk.toString(),sbv.toString());
				skipws();
				if( this.k < this.mapStr.length() && this.mapStr.charAt(this.k)==';')
					{
					this.k++;
					skipws();
					}
				return entry;
				}
			}
		}
 	
	
	}
