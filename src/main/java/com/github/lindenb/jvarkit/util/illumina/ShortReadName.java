/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2015 creation

*/

package com.github.lindenb.jvarkit.util.illumina;

import java.util.Arrays;

import htsjdk.samtools.SAMRecord;

/* http://en.wikipedia.org/wiki/FASTQ_format
 * Parse a Read name. For now only  Casava 1.8 
 */
public class ShortReadName
	{
	private boolean valid=false;
	private String tokens[]=new String[7];
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
		return valid;
		}
	
	public String getInstrumentName()
		{
		return (valid?tokens[0]:null);
		}
	
	public int getRunId()
		{
		return (valid?Integer.parseInt(tokens[1]):-1);
		}
	
	public String getFlowCellId()
		{
		return (valid?tokens[2]:null);
		}
	
	public int getFlowCellLane()
		{
		return (valid?Integer.parseInt(tokens[3]):-1);
		}
	
	public int getTile()
		{
		return (valid?Integer.parseInt(tokens[4]):-1);
		}
	
	public int getX()
		{
		return (valid?Integer.parseInt(tokens[5]):-1);
		}
	public int getY()
		{
		return (valid?Integer.parseInt(tokens[6]):-1);
		}

	@Override
	public int hashCode()
		{
		if(valid)
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
		if (valid != other.valid)
			return false;
		if (!Arrays.equals(tokens, other.tokens))
			return false;
		return true;
		}

	@Override
	public String toString()
		{
		StringBuilder b=new StringBuilder();
		for(int i=0;i< tokens.length;++i)
			{
			if(i>0) b.append(":");
			b.append(tokens[i]);
			}
		return b.toString();
		}
	
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
		
		int col=0;
		int prev=-1;
		while(col < r.tokens.length)
			{
			int colon=readName.indexOf(':', prev+1);
			r.tokens[col] =(colon==-1?
					readName.substring(prev+1):
					readName.substring(prev+1,colon)
					);
			col++;
			if(colon==-1) break;
			prev=colon;
			}
		r.valid=(col==r.tokens.length &&
				positiveInt(r.tokens[1]) &&
				positiveInt(r.tokens[3]) &&
				positiveInt(r.tokens[4]) &&
				positiveInt(r.tokens[5]) &&
				positiveInt(r.tokens[6])
				);
		return r;
		}
	
	/** Parse a Read name. For now only  Casava 1.8 */
	public static ShortReadName parse(SAMRecord record)
		{
		return parse(record==null?null:record.getReadName());
		}
	}
