package com.github.lindenb.jvarkit.util.picard;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import net.sf.picard.PicardException;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

/** in the SAM Spec 2014, XP was renamed SA */ 
public class OtherCanonicalAlignFactory {
	private SAMFileHeader header;
	private final Pattern semiColon=Pattern.compile("[;]");
	private final Pattern comma=Pattern.compile("[,]");
	private String tagname="SA"; 
	public OtherCanonicalAlignFactory(SAMFileHeader header)
		{
		this.header=header;
		}
	
	
	public String getAttributeKey()
		{
		return tagname;
		}
	
	public void setAttributeKey(String tagname)
		{
		if(tagname==null || tagname.length()!=2) throw new IllegalArgumentException("Bad tag name");
		this.tagname = tagname;
		}
	
	public List<OtherCanonicalAlign> getXPAligns(final SAMRecord record)
		{
	
		String xp=record.getStringAttribute(getAttributeKey());
		if(xp==null) return Collections.emptyList();
		
		String ss[]=this.semiColon.split(xp);
		List<OtherCanonicalAlign> L=new ArrayList<OtherCanonicalAlign>(ss.length);
		for(String s:ss)
			{
			if(s.isEmpty()) continue;
			
			String tokens[]=this.comma.split(s);

			XPAlignImpl f2=new XPAlignImpl();
			
			f2.tid= this.header.getSequenceIndex(tokens[0]);
			if(f2.tid==-1) throw new PicardException("Unknown chromosome in "+s);
			f2.pos=Integer.parseInt(tokens[1]);
			f2.strand=tokens[2].charAt(0);
			f2.cigarStr=tokens[3];
			f2.mapq=Integer.parseInt(tokens[4]);	
			f2.nm=Integer.parseInt(tokens[5]);
				
			L.add(f2);
			}
		return L;
		}
	
	
	private class XPAlignImpl  implements OtherCanonicalAlign
		{
		private int tid;
		private int pos;
		private char strand;
		private String cigarStr;
		private Cigar cigar=null;
		private int mapq=0;
		private int nm=0;
		
		@Override
		public int getChromIndex() {
			return tid;
			}
		
		@Override
		public String getChrom()
			{
			return OtherCanonicalAlignFactory.this.header.getSequenceDictionary().getSequence(this.tid).getSequenceName();
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
		public int getMapQ() {
			return mapq;
			}
		
		@Override
		public int getNM() {
			return nm;
			}
		
		@Override
		public int compareTo(OtherCanonicalAlign o) {
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
			result = prime * result + mapq;
			result = prime * result + nm;
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
			OtherCanonicalAlign other = (OtherCanonicalAlign) obj;
			if (pos != other.getPos())
				return false;
			if (strand != other.getStrand())
				return false;
			if (tid != other.getChromIndex())
				return false;
			if (!cigarStr.equals(other.getCigar()))
				return false;
			if(mapq!=other.getMapQ()) return false;
			if(nm!=other.getNM()) return false;
			return true;
			}
	
		@Override
		public String toString() {
			return getChrom()+","+getStrand()+getPos()+","+getCigarString()+","+mapq+","+nm;
			}
		
	
}

	
}
