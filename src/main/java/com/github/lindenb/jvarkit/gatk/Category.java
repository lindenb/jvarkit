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
* 2017 creation

*/
package com.github.lindenb.jvarkit.gatk;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;


public class Category implements Iterable<Object>{
	private final int _hash;
	private final Object labels[] ;
	public Category(final List<Object> labels) {
		this.labels= labels.toArray(new Object[labels.size()]);
		this._hash = Arrays.hashCode(this.labels);
		}
	
	@Override
	public Iterator<Object> iterator() {
		return Arrays.asList(labels).iterator();
		}
	
	public int size() {
		return this.labels.length;
	}
	
	public Object get(int i){
		return this.labels[i];
	}
	
	@Override
	public int hashCode() {
		return this._hash;
		}
	@Override
	public boolean equals(final Object o) {
		if(o==this) return true;
		if(o==null || !(o instanceof Category)) return false;
		final Category other=Category.class.cast(o);
		return this._hash==other._hash &&
				Arrays.equals(this.labels,other.labels);
		}
	
	@Override
	public String toString() {
		return this.labels.toString();
		}
	}
