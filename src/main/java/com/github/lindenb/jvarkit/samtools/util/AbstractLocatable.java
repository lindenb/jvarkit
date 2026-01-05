/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.samtools.util;


import htsjdk.samtools.util.Locatable;

public abstract class AbstractLocatable implements ExtendedLocatable {

@Override
public boolean equals(final Object obj) {
	if(this==obj) return true;
	if(obj==null || !(obj instanceof AbstractLocatable)) return false;
	final AbstractLocatable o = AbstractLocatable.class.cast(obj);
	if(this.getStart()!=o.getStart()) return false;
	if(this.getEnd()!=o.getEnd()) return false;
	return this.contigsMatch(o);
	}

@Deprecated
public static int compareTo(final Locatable A, Locatable B) {
	return LocatableUtils.compareTo(A, B);
	}

@Override
public int hashCode() {
	final int prime = 31;
	int result = 1;
	result = prime * result + getContig().hashCode();
	result = prime * result + Integer.hashCode(getStart());
	result = prime * result + Integer.hashCode(getEnd());
	return result;	
}

@Override
public String toString() {
	return getContig()+":"+getStart()+"-"+getEnd();
	}
}
