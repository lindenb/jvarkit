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
import java.util.Collections;
import java.util.Arrays;

public interface Sample {
	public String getId();
	public Family getFamily();
	public Status getStatus();
	public Sex getSex();
	
	public default boolean isMale() { return Sex.male.equals(this.getSex());}
	public default boolean isFemale() { return Sex.female.equals(this.getSex());}

	public default boolean isAffected() { return Status.affected.equals(this.getStatus());}
	public default boolean isUnaffected() { return Status.unaffected.equals(this.getStatus());}
        
	public boolean hasFather();
	public boolean hasMother();
	public Sample getFather();
	public Sample getMother();
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
	}
