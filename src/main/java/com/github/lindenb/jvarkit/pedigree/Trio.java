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
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/** a trio children/father/mother */
public interface Trio extends SampleSet, Iterable<Sample>, Comparable<Trio> {
	public default boolean hasFather() { return getFather()!=null;}
	public default boolean hasMother() { return getMother()!=null;}
	public default Sample getFather() { return getChild().getFather();}
	public default Sample getMother() { return getChild().getMother();}
	
	/** return the parents of the child */
	public default List<Sample> getParents() {
		return getChild().getParents();
	}
	
	/** child of this trio */
	public Sample getChild();
	
	@Override
	public default Iterator<Sample> iterator() {
		return getSamples().iterator();
		}
	
	@Override
	public default Set<Sample> getSamples() {
		final Set<Sample> L = new TreeSet<>();
		if(hasFather()) L.add(getFather());
		if(hasMother()) L.add(getMother());
		L.add(getChild());
		return L;
		}
	@Override
	public default int compareTo(final Trio o) {
		return getChild().compareTo(o.getChild());
		}
	}
