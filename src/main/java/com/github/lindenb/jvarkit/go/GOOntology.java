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
package com.github.lindenb.jvarkit.go;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;



/**
 * Gene Ontology tree
 * @author lindenb
 *
 */
public interface GOOntology extends Iterable<GOOntology.Term>
	{

	public static final String GO_OBO_URL="http://current.geneontology.org/ontology/go-basic.obo";
	
	public static interface Term {
		public Set<String> getSynonyms();
		public Division getDivision();
		public String getAcn();
		public String getName();
		public String getDefinition();
		public boolean hasRelations();
		public Set<Relation> getRelations();
		public boolean isDescendantOf(final Term parentNode);
		public int getMinDepth();
	}

	public static enum Division
		{
		cellular_component("GO:0005575"),
		biological_process("GO:0008150"),
		molecular_function("GO:0003674");
		private final String acn;
		Division(final String acn) { this.acn=acn;}
		/** get accession number */
		public String getAcn() { return this.acn;}
		}

	
	public static interface Relation
		{
		public String getType();
		public Term getTo();
		}
	
	
	
	
	
	public int size();
	public Collection<GOOntology.Term> getTerms();
	@Override
	public Iterator<GOOntology.Term> iterator();
	

	
	public Term getTermByAccession(final String s);
	
	/** search term by name or accession or synonyms, ignoring case */
	public Term getTermByName(final String s);
	
	public default Term getTermByAccessionOrName(final String s) {
		Term t= this.getTermByAccession(s);
		if(t==null)
			{
			t= this.getTermByName(s);
			}
		return t;
		}

	}
