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
package com.github.lindenb.jvarkit.lang;
import java.math.BigInteger;
import java.util.Comparator;

/**
 * SmartComparator. Default is NOT case sensitive
 *
 */
public class SmartComparator implements Comparator<String> {
	private boolean case_sensitive = false;
	
	public SmartComparator() {
		}
	private SmartComparator(final SmartComparator cp) {
		this.case_sensitive = cp.case_sensitive;
	}
	
	/** clone the comparator and return a case-sensitive one */
	public SmartComparator caseSensitive() {
		final SmartComparator copy = new SmartComparator(this);
		copy.case_sensitive = true;
		return copy;
		}
	
	public boolean isCaseSensitive() {
		return case_sensitive;
		}
	
	private char tr(final char c) {
		return this.case_sensitive ? c : Character.toUpperCase(c);
		}
	
	@Override
    public int compare(final String a, final String b)
        {
        if(a==null && b==null) return 0;
        if(a==null && b!=null) return -1;
        if(a!=null && b==null) return 1;
        final String ss[]={a,b};
        final char c[]={0,0};
        final int i[]={0,0};
        for(;;)
            {
            if(i[0] == a.length() &&
               i[1] == b.length()) return 0;
            if(i[0] == a.length()) return -1;
            if(i[1] == b.length()) return 1;

            c[0] = tr(a.charAt(i[0]));
            c[1] = tr(b.charAt(i[1]));

            //search for integer
            if( Character.isDigit(c[0]) &&
                Character.isDigit(c[1]))
                {
                int j[]={i[0]+1,i[1]+1};
                try
	                {
	                final BigInteger values[]={null,null};
	                for(int side=0;side<2;++side)
	                    {
	                    while(  j[side]< ss[side].length() &&
	                                Character.isDigit(ss[side].charAt(j[side])))
                            {
                            j[side]++;
                            }
	                    values[side]=new BigInteger(ss[side].substring(i[side],j[side]));
	                    }
	                final int k = values[0].compareTo(values[1]);
	                i[0]=j[0];
	                i[1]=j[1];
	                if(k!=0) return k;
	                continue;
	                }
                catch(final NumberFormatException err)
                	{
                	throw err;
                	}
                }           
            int k=c[0]-c[1];
            if(k!=0) return (k<0?-1:1);
            i[0]++;
            i[1]++;
            }
        }
	}
