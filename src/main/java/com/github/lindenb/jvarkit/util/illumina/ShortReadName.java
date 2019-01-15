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

package com.github.lindenb.jvarkit.util.illumina;

import java.util.Arrays;
import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.samtools.SAMRecord;

/* http://en.wikipedia.org/wiki/FASTQ_format
 * Parse a Read name. For now only  Casava 1.8 
 */
public class ShortReadName
	{
	private enum Type {INVALID,SEVEN_COLS,SIX_COLS};
	private Type type = Type.INVALID;
	private String tokens[];
	private ShortReadName()
		{
		
		}
	
	private static boolean positiveInt(String s)
		{
		try {
			Integer.parseInt(s);
			return true;
			} 
		catch (NumberFormatException e) {
			return false;
			}
		}
	
	
	public boolean isValid()
		{
		return !this.type.equals(Type.INVALID);
		}
	
	public String getInstrumentName()
		{
		switch(type)
			{
			case INVALID : return null;
			case SEVEN_COLS : return tokens[0];
			case SIX_COLS: return tokens[0];
			default: throw new IllegalStateException();
			}
		}
	
	public int getRunId()
		{
		switch(type)
			{
			case INVALID : return -1;
			case SEVEN_COLS : return Integer.parseInt(tokens[1]);
			case SIX_COLS : return -1;
			default: throw new IllegalStateException();
			}
		}
	
	public String getFlowCellId()
		{
		switch(type)
			{
			case INVALID : return null;
			case SEVEN_COLS : return tokens[2];
			case SIX_COLS : return tokens[1];
			default: throw new IllegalStateException();
			}
		}
	
	public int getFlowCellLane()
		{
		switch(type)
			{
			case INVALID : return -1;
			case SEVEN_COLS : return Integer.parseInt(tokens[3]);
			case SIX_COLS :return  Integer.parseInt(tokens[2]);
			default: throw new IllegalStateException();
			}
		}
	
	public int getTile()
		{
		switch(type)
			{
			case INVALID : return -1;
			case SEVEN_COLS : return Integer.parseInt(tokens[4]);
			case SIX_COLS :  return Integer.parseInt(tokens[3]);
			default: throw new IllegalStateException();
			}
		}
	
	public int getX()
		{
		switch(type)
			{
			case INVALID : return -1;
			case SEVEN_COLS : return Integer.parseInt(tokens[5]);
			case SIX_COLS : return Integer.parseInt(tokens[4]);
			default: throw new IllegalStateException();
			}
		}
	public int getY()
		{
		switch(type)
			{
			case INVALID : return -1;
			case SEVEN_COLS : return Integer.parseInt(tokens[6]);
			case SIX_COLS : return Integer.parseInt(tokens[5]);
			default: throw new IllegalStateException();
			}
		}

	@Override
	public int hashCode()
		{
		if(isValid())
			{
			return 31 + Arrays.hashCode(tokens);
			}
		else
			{
			return 31;
			}
		}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ShortReadName other = (ShortReadName) obj;
		if (isValid() != other.isValid())
			return false;
		if (!Arrays.equals(tokens, other.tokens))
			return false;
		return true;
		}
	
	public String getName()
		{
		return String.join(":", this.tokens);
		}

	
	@Override
	public String toString()
		{
		return this.getName();
		}
	
	private static final CharSplitter COLON=CharSplitter.COLON;
	
	/** Parse a Read name. For now only  Casava 1.8 */
	public static ShortReadName parse(String readName)
		{
		
		ShortReadName r=new ShortReadName();
		
		if(readName==null || readName.isEmpty()) return r;
		if(readName.charAt(0)=='@') readName=readName.substring(1);
		if(readName.endsWith("/1") || readName.endsWith("/2"))
			readName=readName.substring(0,readName.length()-2);
		int ws= readName.indexOf(' ');
		if(ws!=-1) ws= readName.indexOf('\t');
		if(ws!=-1)  readName=readName.substring(0,ws);
		
		r.tokens = COLON.split(readName);
		if(	r.tokens.length == 7  &&
				positiveInt(r.tokens[1]) &&
				positiveInt(r.tokens[3]) &&
				positiveInt(r.tokens[4]) &&
				positiveInt(r.tokens[5]) &&
				positiveInt(r.tokens[6])
				)
			{
			r.type = Type.SEVEN_COLS;
			}
		else if(	r.tokens.length == 6  &&
				positiveInt(r.tokens[2]) &&
				positiveInt(r.tokens[3]) &&
				positiveInt(r.tokens[4]) &&
				positiveInt(r.tokens[5])
				)
			{
			r.type = Type.SIX_COLS;
			}
		return r;
		}
	
	/** Parse a Read name. For now only  Casava 1.8 */
	public static ShortReadName parse(final SAMRecord record)
		{
		return parse(record==null?null:record.getReadName());
		}
	}
