package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

/**  utilities for BWA */
public class BWAUtils
	{
	public static class XPAlign implements Comparable<XPAlign>
		{
		private String chrom;
		private int pos;
		private char strand;
		private String cigarStr;
		private Cigar cigar=null;
		public String getChrom()
			{
			return chrom;
			}
		public int getPos()
			{
			return pos;
			}
		public char getStrand()
			{
			return strand;
			}
		public String getCigarString()
			{
			return cigarStr;
			}
		public Cigar getCigar()
			{
			if(cigar==null) cigar=TextCigarCodec.getSingleton().decode(getCigarString());
			return cigar;
			}
		public List<CigarElement> getCigarElements()
			{
			return getCigar().getCigarElements();
			}
		
		@Override
		public int compareTo(XPAlign o) {
			int i=getChrom().compareTo(o.getChrom());
			if(i!=0) return i;
			i=getPos()-(o.getPos());
			if(i!=0) return i;
			i=getStrand()-(o.getStrand());
			if(i!=0) return i;
			return getCigarString().compareTo(o.getCigarString());
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result +chrom.hashCode();
			result = prime * result + cigar.hashCode();
			result = prime * result + cigarStr.hashCode();
			result = prime * result + pos;
			result = prime * result + strand;
			return result;
			}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			XPAlign other = (XPAlign) obj;
			if (pos != other.pos)
				return false;
			if (strand != other.strand)
				return false;
			if (!chrom.equals(other.chrom))
				return false;
			if (!cigarStr.equals(other.cigarStr))
				return false;
			return true;
			}

		@Override
		public String toString() {
			return getChrom()+","+getStrand()+getPos()+","+getCigarString();
			}
		}
	
	public static List<XPAlign> getXPAligns(SAMRecord record)
		{
	
		String xp=record.getStringAttribute("XP");
		if(xp==null) return Collections.emptyList();
		List<XPAlign> L=new ArrayList<XPAlign>();
		
		Pattern semiColon=Pattern.compile("[;]");
		
		Pattern comma=Pattern.compile("[,]");
		for(String s:semiColon.split(xp))
			{
			if(s.isEmpty()) continue;
			
			String tokens[]=comma.split(s);

			XPAlign f2=new XPAlign();
			f2.chrom=tokens[0];
			f2.pos=Integer.parseInt(tokens[1].substring(1));
			f2.strand=tokens[1].charAt(0);
			f2.cigarStr=tokens[2];
	
			L.add(f2);
			}
		return L;
		}

}
