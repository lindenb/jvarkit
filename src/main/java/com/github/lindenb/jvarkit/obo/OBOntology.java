/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.obo;

import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/** direct graph for Open Biological and Biomedical Ontology  (OBO) */
public interface OBOntology extends Iterable<OBOntology.Term> {

	/** a term in the ontology */
	public static interface Term
		{
		public String getAcn();
		public String getLabel();
		/** get description optionaly read by the OBO parser */
		public Optional<String> getDescription();
		/* get direct parents */
		public Set<Term> getParents();
		/* get direct children */
		public Set<Term> getChildren();
		/* get all ancestors */
		public Set<Term> getAllAncestors();
		/* get direct children */
		public Set<Term> getAllDescendant();
		/* return true if this==node or this is a parent of node */
		public boolean isAncestorOf(final Term node);
		/* return true if this==node or this is a parent of node */
		public boolean isDescendantOf(final Term node);
		}
	/** find a term by Acn or return null */
	public  Term findTermByAcn(final String id);
	
	/** find a term by Acn or name, or return null */
	public  default Term findTermByAcnOrName(final String id) {
		final Term t = findTermByAcn(id);
		if(t!=null) return t;
		final String s1 = id.replace(' ', '_');
		final String s2 = id.replace('_', ' ');
		return stream().
				filter(T->T.getLabel().equalsIgnoreCase(s1) || T.getLabel().equals(s2)).
				findAny().orElse(null);
		}

	
	/** return a stream of the terms of this ontology */
	public Stream<Term> stream();
	
	/** return the stream of Terms as List */
	public default List<Term> toList() {
		return stream().collect(Collectors.toList());
	}
	
	@Override
	public default Iterator<Term> iterator() {
		return stream().iterator();
		}
	}
