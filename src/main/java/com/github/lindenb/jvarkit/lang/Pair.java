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
package com.github.lindenb.jvarkit.lang;

import java.util.Map;

public class Pair<T1,T2> {
	protected final T1 first;
	protected final T2 second;
	public Pair(final Map.Entry<T1, T2> p) {
		this(p.getKey(),p.getValue());
		}
	
	public Pair(final T1 first,final T2 second) {
		this.first = first;
		this.second = second;
		}

	protected boolean sameAny(final Object a, final Object b) {
		if(a==b) return true;
		if(a==null && b==null) return true;
		if(a==null || b==null) return false;
		if(!a.getClass().equals(b.getClass())) return false;
		return a.equals(b);
		}
	
	protected boolean sameFirst(final Object a, final Object b) {
		return sameAny(a,b);
		}
	protected boolean sameSecond(final Object a, final Object b) {
		return sameAny(a,b);
		}
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof Pair)) return false;
		final Pair<?,?> p=(Pair<?,?>)obj;
		return sameFirst(this.getFirst(),p.getFirst()) &&
				sameSecond(this.getSecond(), p.getSecond());
		}
	public T1 getFirst() {
		return first;
		}
	public T2 getSecond() {
		return second;
		}
	@Override
	public int hashCode() {
		return (first==null?-1:first.hashCode())*31 +
					(second==null?-1:second.hashCode());
		}
	@Override
	public String toString() {
		return "("+getFirst()+","+getSecond()+")";
		}
	}
