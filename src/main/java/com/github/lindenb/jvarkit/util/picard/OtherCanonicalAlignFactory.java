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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

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
			if(f2.tid==-1) throw new RuntimeException("Unknown chromosome in "+s);
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
		private int  mAlignmentEnd = SAMRecord.NO_ALIGNMENT_START;
		
		@Override
		public int getChromIndex() {
			return tid;
			}
		
		@Override
		public String getReferenceName() {
			return OtherCanonicalAlignFactory.this.header.getSequenceDictionary().getSequence(this.tid).getSequenceName();
			}	
		
		@Override
		public String getChrom()
			{
			return getReferenceName();
			}
		
		@Override
		public int getAlignmentStart() {
			return pos;
			}
		@Override
		public int getPos()
			{
			return getAlignmentStart();
			}
		@Override
		public boolean getReadNegativeStrandFlag() {
			return strand=='-';
			}
		@Override
		public char getStrand()
			{
			return getReadNegativeStrandFlag()?'-':'+';
			}

		@Override
		public String getCigarString()
			{
			return cigarStr;
			}
		@Override
		public Cigar getCigar()
			{
			if(cigar==null) cigar= TextCigarCodec.decode(getCigarString());
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
			i=getAlignmentStart()-(o.getAlignmentStart());
			if(i!=0) return i;
			i=(getReadNegativeStrandFlag()?-1:1)-(o.getReadNegativeStrandFlag()?-1:1);
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
			if (pos != other.getAlignmentStart())
				return false;
			if ((strand=='-') != other.getReadNegativeStrandFlag())
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
			return getReferenceName()+","+getStrand()+getAlignmentStart()+","+getCigarString()+","+mapq+","+nm;
			}
		
		@Override
	    public int getUnclippedStart() {
	        int pos = this.getPos();

	        for (final CigarElement cig : getCigar().getCigarElements()) {
	            final CigarOperator op = cig.getOperator();
	            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
	                pos -= cig.getLength();
	            }
	            else {
	                break;
	            }
	        }

	        return pos;
	    	}

		@Override
		public int getUnclippedEnd() {
		        int pos = getAlignmentEnd();
		        final List<CigarElement> cigs = getCigar().getCigarElements();
		        for (int i=cigs.size() - 1; i>=0; --i) {
		            final CigarElement cig = cigs.get(i);
		            final CigarOperator op = cig.getOperator();

		            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
		                pos += cig.getLength();
		            }
		            else {
		                break;
		            }
		        }

		        return pos;               
		    }
		

	    /**
	     * @return 1-based inclusive rightmost position of the clipped sequence, or 0 read if unmapped.
	     */
		@Override
	    public int getAlignmentEnd() {
	       if(this.mAlignmentEnd == SAMRecord.NO_ALIGNMENT_START) {
	            this.mAlignmentEnd = getPos() + getCigar().getReferenceLength() - 1;
	       		}
	        return this.mAlignmentEnd;
	    	}
		
	
		}

	
}
