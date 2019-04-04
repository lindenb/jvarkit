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
package com.github.lindenb.jvarkit.pedigree;
import java.util.Collections;
import java.util.Set;

class TrioImpl implements Trio {
	private final Sample child;
	TrioImpl(final Sample child)
		{
		this.child=child;
		if(child==null) throw new IllegalArgumentException("child==null");
		if(!child.hasFather() && !child.hasMother() ) throw new IllegalArgumentException("child without parents "+child);
		}
	@Override
	public Sample getChild() { return this.child;}
	
	@Override
	public Set<Trio> getTrios() {
		return Collections.singleton(this);		
		}
	@Override
	public int hashCode() {
		return getChild().hashCode();	
		}
	@Override
	public boolean equals(final Object o)
		{
		if(o==this) return true;
		if(o==null || !(o instanceof TrioImpl)) return false;
		return getChild().equals(TrioImpl.class.cast(o).getChild());
		}
	@Override
	public String toString() {
		return "trio(child:"+getChild()+",father:"+ getFather()+",mother:"+getMother()+")";		
		}
	}
