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
