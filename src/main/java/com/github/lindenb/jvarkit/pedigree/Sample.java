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
import java.util.List;

import com.github.lindenb.jvarkit.lang.StringUtils;

import java.util.Collections;
import java.util.Arrays;

/** a sample in a family in a pedigree */
public interface Sample extends Comparable<Sample> {
	public String getId();
	public Family getFamily();
	public Status getStatus();
	public Sex getSex();
	
	
	public default Pedigree getPedigree() { return getFamily().getPedigree();}

	
	public default boolean isMale() { return Sex.male.equals(this.getSex());}
	public default boolean isFemale() { return Sex.female.equals(this.getSex());}

	public default boolean isAffected() { return Status.affected.equals(this.getStatus());}
	public default boolean isUnaffected() { return Status.unaffected.equals(this.getStatus());}
        
	public default boolean hasFather() { return getFather()!=null;}
	public default boolean hasMother() { return getMother()!=null;}
	/** return the father of this sample, or null */
	public Sample getFather();
	/** return the mother of this sample, or null */
	public Sample getMother();
	
	/** get parent using index 0: father, 1: mother, other values will throw an exception */
	public default Sample getParent( int zeroOrOne) {
		switch(zeroOrOne) {
		case 0: return getFather();
		case 1: return getMother();
		default: throw new IllegalArgumentException("0 or 1 but got "+zeroOrOne);
		} 
	}
	
	
	
	/** get the parents of this individual */
	public default List<Sample> getParents() {
		final Sample f = getFather();
		final Sample m = getFather();
		if(f==null)
			{
			if(m==null) return Collections.emptyList();
			return Collections.singletonList(m);
			}
		else {
			if(m!=null) return Arrays.asList(f,m);
			return Collections.singletonList(f);
			}
		}
	
	/** return true if sample has unique id in the pedigree */
	public default boolean hasUniqId() 
		{
		return getPedigree().
				getSamples().
				stream().
				filter(P->getId().equals(P.getId())).
				limit(2L).
				count() ==1L;
		}
	
	@Override
	default int compareTo(final Sample o) {
		int i = getFamily().compareTo(o.getFamily());
		if(i!=0) return i;
		return getId().compareTo(o.getId());
		}
	}
