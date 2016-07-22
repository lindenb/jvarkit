/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.picard;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalUtil;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class IntervalUtils extends IntervalUtil
	{
	public static Interval parseOne(String s)
		{
		return parseOne(null,s);
		}
	
	public static Interval parseOne(SAMSequenceDictionary dict,String s)
		{
		int colon=s.indexOf(':');
		if(colon<1 || colon+1==s.length()) return null;
		int hyphen=s.indexOf('-',colon+1);
		if(hyphen==-1) return null;
		int start,end;
		try {
			SAMSequenceRecord ssr=null;
			String chrom=s.substring(0,colon);
			if(dict!=null)
				{
				ssr=dict.getSequence(chrom);
				if(ssr==null) return null;
				}
			start=Integer.parseInt(s.substring(colon+1,hyphen));
			end=Integer.parseInt(s.substring(hyphen+1));
			if(end<start) return null;
			return new Interval(s.substring(0,colon), start, end);
		} catch (Exception e)
			{
			return null;
			}
		
		}
	}
