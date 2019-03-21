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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.lang;

/** An abstract char sequence */
public abstract class AbstractCharSequence
	implements CharSequence
	{
	@Override
	public int hashCode()
		{
        int hash=0;
        int L=this.length();
        for (int i = 0; i < L; i++) {
        	hash = 31 * hash + (int)charAt(i);
            }
        return hash;
		}
	
	public boolean hasStringValue(final String s)
		{
		if(s.length()!=this.length()) return false;
		for(int i=0;i< s.length();++i)
			{
			if(this.charAt(i)!=s.charAt(i)) return false;
			}
		return true;
		}
	
	public String getString()
		{
		final StringBuilder b=new StringBuilder(length());
		for(int i=0;i< length();++i) b.append(charAt(i));
		return b.toString();
		}
	
	@Override
	public String toString()
		{
		return getString();
		}
	@Override
	public CharSequence subSequence(final int start,final int end)
		{
		if(start==0 && end==this.length()) return this;
		return new SubSequence(this,start,end);
		}
	}
