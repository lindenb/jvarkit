package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

public class XPalignFactory {
	private SAMFileHeader header;
	private final Pattern semiColon=Pattern.compile("[;]");
	private final Pattern comma=Pattern.compile("[,]");

	public XPalignFactory(SAMFileHeader header)
		{
		this.header=header;
		}
	
	
	public List<XPAlign> getXPAligns(final SAMRecord record)
		{
	
		String xp=record.getStringAttribute("XP");
		if(xp==null) return Collections.emptyList();
		
		String ss[]=this.semiColon.split(xp);
		List<XPAlign> L=new ArrayList<XPAlign>(ss.length);
		for(String s:ss)
			{
			if(s.isEmpty()) continue;
			
			String tokens[]=this.comma.split(s);
	
			XPAlignImpl f2=new XPAlignImpl();
			f2.tid=record.getReferenceIndex();
			f2.pos=Integer.parseInt(tokens[1].substring(1));
			f2.strand=tokens[1].charAt(0);
			f2.cigarStr=tokens[2];
			L.add(f2);
			}
		return L;
		}
	
	
	private class XPAlignImpl  implements XPAlign
		{
		private int tid;
		private int pos;
		private char strand;
		private String cigarStr;
		private Cigar cigar=null;
		
		@Override
		public int getChromIndex() {
			return tid;
			}
		
		@Override
		public String getChrom()
			{
			return XPalignFactory.this.header.getSequenceDictionary().getSequence(this.tid).getSequenceName();
			}
		@Override
		public int getPos()
			{
			return pos;
			}
		@Override
		public char getStrand()
			{
			return strand;
			}
		@Override
		public String getCigarString()
			{
			return cigarStr;
			}
		@Override
		public Cigar getCigar()
			{
			if(cigar==null) cigar=TextCigarCodec.getSingleton().decode(getCigarString());
			return cigar;
			}
		@Override
		public List<CigarElement> getCigarElements()
			{
			return getCigar().getCigarElements();
			}
		
		@Override
		public int compareTo(XPAlign o) {
			int i= this.tid - o.getChromIndex();
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
			result = prime * result +tid;
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
			if (pos != other.getPos())
				return false;
			if (strand != other.getStrand())
				return false;
			if (tid != other.getChromIndex())
				return false;
			if (!cigarStr.equals(other.getCigar()))
				return false;
			return true;
			}
	
		@Override
		public String toString() {
			return getChrom()+","+getStrand()+getPos()+","+getCigarString();
			}
		
	
}

	
}
