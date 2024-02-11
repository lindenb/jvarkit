package com.github.lindenb.jvarkit.samtools.util;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;

public class LocatableUtils extends CoordMath {
	public static int compareTo(final Locatable A, Locatable B) {
		int i= A.getContig().compareTo(B.getContig());
		if(i!=0) return i;
		i = Integer.compare(A.getStart(),B.getStart());
		if(i!=0) return i;
		i = Integer.compare(A.getEnd(),B.getEnd());
		return i;
		}
	
	/** use StringUtils.niceInt to format start and end  */
	public static String toNiceString(final Locatable loc) {
		return loc.getContig()+":"+ 
				StringUtils.niceInt(loc.getStart()) + "-"+ 
				StringUtils.niceInt(loc.getEnd())
				;
		}
	/** use StringUtils.niceInt to format start and end  */
	public static String toBed3(final Locatable loc) {
		return String.join(
				"\t",
				loc.getContig(),
				String.valueOf(loc.getStart()-1) ,
				String.valueOf(loc.getEnd())
				);
		}
}
