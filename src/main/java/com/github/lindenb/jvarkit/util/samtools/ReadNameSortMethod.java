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
package com.github.lindenb.jvarkit.util.samtools;

import java.util.Comparator;
import java.util.function.Supplier;

import htsjdk.samtools.SAMRecord;

/** comparator for sam sorted on query-name (not the same between picard and samtools sort */
public enum ReadNameSortMethod implements Supplier<Comparator<SAMRecord>> {
samtools
	{
	@Override
	public Comparator<SAMRecord> get() {
		return SAMTOOLS_CMP;
		}
	},
picard{
	@Override
	public Comparator<SAMRecord> get() {
		return PICARD_CMP;
		}
	}
;
	
public static final String DESCRIPTION="Method used to sort the read on query name. (samtools != picard) see https://github.com/samtools/hts-specs/issues/5";


private static int side(final SAMRecord rec) {
	if(rec.getReadPairedFlag()) {
		if(rec.getFirstOfPairFlag()) return 1;
		if(rec.getSecondOfPairFlag()) return 2;
		throw new IllegalArgumentException("Side for flag "+rec.getReadName()+":"+rec.getFlags()+"?");
	}
	else {
		return 0;
	}
}

private static  final Comparator<SAMRecord> SAMTOOLS_CMP = new SamToolsReadNameComparator();
private static final  Comparator<SAMRecord> PICARD_CMP = new SimpleReadNameComparator();

/** see https://github.com/samtools/hts-specs/issues/5 */
private static class SamToolsReadNameComparator implements Comparator<SAMRecord>
	{
	@Override
	public int compare(final SAMRecord o1, final SAMRecord o2) {
		int i= strnum_cmp(o1.getReadName(), o2.getReadName());
		if(i!=0) return i;
		return side(o1)-side(o2);

		}
	
	/* see https://github.com/samtools/hts-specs/issues/5 */
	private int strnum_cmp(final String _a,final String _b)
		{
	
			char ca='\0',cb='\0';
		    int ia = 0, ib = 0;
		    while(ia< _a.length() && ib<_b.length()) {
		    	ca=(ia< _a.length()?_a.charAt(ia):'\0');
		    	cb=(ib< _b.length()?_b.charAt(ib):'\0');
		    	
		    	
		        if (Character.isDigit(ca) && Character.isDigit(cb))
		        	{
		            while (ca=='0'){
		            	++ia;
		            	ca=(ia< _a.length()?_a.charAt(ia):'\0');
		            	}
		            while (cb=='0'){
		            	++ib;
		            	cb=(ib< _b.length()?_b.charAt(ib):'\0');
		            	}
		            
		            
		            while (Character.isDigit(ca) && Character.isDigit(cb) && ca==cb)
		            	{
		            	++ia;
		            	++ib;
				    	ca=(ia< _a.length()?_a.charAt(ia):'\0');
				    	cb=(ib< _b.length()?_b.charAt(ib):'\0');
		            	}
		            
		            if (Character.isDigit(ca) && Character.isDigit(cb)) {
		                int i = 0;
		                while((ia+i)< _a.length() &&
		                	  (ib+i)< _b.length() &&
		                	  Character.isDigit(_a.charAt(ia+i)) &&
		                	  Character.isDigit(_b.charAt(ib+i))
		                	  ) {
		                	++i;
		                	}
				    	final char c1=((ia+i)< _a.length()?_a.charAt(ia+i):'\0');
				    	final char c2=((ib+i)< _b.length()?_b.charAt(ib+i):'\0');
		                return Character.isDigit(c1)? 1 : Character.isDigit(c2)? -1 : (int)ca - (int)cb;
		            }
		          else if (Character.isDigit(ca))
		        	  {
		        	  return 1;
		        	  }
		         else if (Character.isDigit(cb))
		        	 {
		        	 return -1;
		        	 }
		         else if (ia != ib)
		        	 {
		        	 return ia - ib;
		        	 }
		        	}/* end of is digit */
		        else
		        	{
		            if (ca != cb) {
		            	return (int)ca - (int)cb;
		            }
		            ++ia; ++ib;
		        	}
		    	}
		    if(ca==cb) return 0;
		    return ca=='\0'?1:cb=='\0'?-1: (int)ca - (int)cb;
		}
	}

private static class SimpleReadNameComparator implements Comparator<SAMRecord>
	{
	@Override
	public int compare(final SAMRecord o1, final SAMRecord o2) {
		int i=o1.getReadName().compareTo(o2.getReadName());
		if(i!=0) return i;
		return side(o1)-side(o2);
		}
}

}
