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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.go.GOOntology.Term;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.StringUtil;

public class GOParser {
	private static Logger LOG=Logger.of(GOParser.class); 
	public static final String GO_URL_OPT_DESC="Gene ontology URI. Formatted as OBO format. ";

	private boolean ignore_synonyms=false;
	private boolean ignore_definitions=false;
	//private boolean ignore_dbxref=false;
	private boolean debug=false;
	
	private static class OBORecord {
		final List<String> lines = new ArrayList<>();
		
		private boolean isHeader(final String s) {
			return s.startsWith("[") && s.endsWith("]");
		}
		private String removeComment(String s) {
			int i= s.indexOf("!");
			if(i==-1) return s.trim();
			return s.substring(0, i).trim();
		}
		private String unquote(String s) {
			if(s.startsWith("\"")) {
				int j = s.lastIndexOf('"');
				if(j!=-1) s=s.substring(1, j);
				}
			return s;
		}
			
		private Map.Entry<String, String> split(String s) {
			int i= s.indexOf(":");
			if(i<=0) throw new IllegalArgumentException("no colon in "+s);
			return new AbstractMap.SimpleEntry<>(
					s.substring(0, i).trim(),
					s.substring(i+1).trim()
					);
		}
		
		List<String> getAll(final String k) {
			return lines.stream().skip(1L).
					map(S->split(S)).
					filter(K->K.getKey().equals(k)).
					map(K->K.getValue()).
					collect(Collectors.toList());
		}
		
		/* get required */
		String get(String k) {
			final List<String> L = getAll(k);
			if(L.size()==1) return L.get(0);
			throw new IllegalArgumentException("Expected one and only one value for "+k+" but got "+L);
		}
		
		/* get optional */
		Optional<String> optional(String k) {
			final List<String> L = getAll(k);
			if(L.size()==1) return Optional.of(L.get(0));
			if(L.isEmpty()) return Optional.empty();
			throw new IllegalArgumentException("Expected one or zero value for "+k+" but got "+L);
		}
		
		boolean has(final String k) {
			return lines.stream().skip(1L).
					map(S->split(S)).
					anyMatch(K->K.getKey().equals(k))
					;
			}
		
		boolean isObsolete() {
			Optional<String> opt = optional("is_obsolete");
			return opt.isPresent() && opt.get().equals("true");
		}
		
		public void add(final String line)  {
			if(lines.isEmpty()) {
				if(!isHeader(line)) throw new IllegalArgumentException("Expected [xxxx] but got "+line);
			} else {
				final String key = split(line).getKey();
				if(key.equals("comment")) return;
				if(key.equals("subset")) return;
				if(key.equals("xref")) return;
				if(key.equals("xref_analog")) return;
				if(key.equals("xref_unk")) return;
				if(key.equals("alt_id")) return;
				if(key.equals("exact_synonym")) return;
				if(key.equals("narrow_synonym")) return;
				if(key.equals("broad_synonym")) return;
				if(key.equals("user_term")) return;
				if(key.equals("created_by")) return;
				if(key.equals("creation_date")) return;
			}
			this.lines.add(line);
		}
		String getId() { return get("id"); }
		String getDefinition() { return unquote(get("def")); }
		String getName() { return unquote(get("name")); }
		String getNamespace() { return get("namespace"); }
		GOOntology.Division getDivision() {
			try {
				return GOOntology.Division.valueOf(getNamespace());
			} catch(final Throwable err) {
				throw new IllegalArgumentException("Cannot convert ["+getNamespace()+"]", err);
			}
		}
		Set<String> getSynonyms() { return getAll("synonym").stream().map(S->unquote(S)).collect(Collectors.toSet()); }
		Set<String> getParents() { return getAll("is_a").stream().map(S->removeComment(S)).collect(Collectors.toSet());}
		boolean isValid() {
			return this.lines.size()>1 && this.lines.get(0).equals("[Term]") && !isObsolete() && has("id");
		}
	}
	
	public GOParser() {}
	
	
	public GOParser setIgnoreDefinitions(boolean ignore_definitions) {
		this.ignore_definitions = ignore_definitions;
		return this;
	}
	
	public GOParser setIgnoreSynonyms(boolean ignore_synonyms) {
		this.ignore_synonyms = ignore_synonyms;
		return this;
		}
	
	public GOParser setDebug(boolean d) {
		this.debug = d;
		return this;
		}
	
	public GOOntology parseOBO(final BufferedReader br) throws IOException
		{
		final GOOntologyImpl tree=new GOOntologyImpl();
		OBORecord rec = new OBORecord();
		boolean in_header= true;
		while(true)
			{
			final String line = br.readLine();
			if(line==null || StringUtil.isBlank(line)) {
				if(rec.isValid()) {
					final TermImpl term = tree.getOrCreate(rec.getId());
					if(!this.ignore_definitions) {
						term.definition = rec.getDefinition();
					}
					term.name = rec.getName();
					term.division = rec.getDivision();
					final Set<String> syns = rec.getSynonyms();
					if(syns!=null && !syns.isEmpty() && !this.ignore_synonyms) {
						term.synonyms = syns;
						}
					for(final String is_a: rec.getParents()) {
						if(is_a.equals(term.getAcn())) continue;
						term.relations.add(new RelationImpl(term, "is_a", tree.getOrCreate(is_a)));
						}
					for(String s: rec.getAll("relationship")) {
						s=rec.removeComment(s);
						final String tokens[] = s.split("\\s+");
						if(tokens.length!=2) throw new IllegalArgumentException("expected two String in "+s);
						if(tokens[1].equals(term.getAcn())) continue;
						term.relations.add(new RelationImpl(term, tokens[0], tree.getOrCreate(tokens[1])));
						}
					}
				if(line==null) break;
				rec = new OBORecord();
				continue;
				}
			
			if(line.equals("[Term]")) in_header=false;
			if(in_header) continue;
			rec.add(line);
			}
						
		
		if(this.debug)
			{
			LOG.debug("tree size: "+tree.size());
			}
		return tree;
		}
	
	
	public GOOntology parseOBO(final File file) throws IOException {
		return parseOBO(file.toPath());
	}

	public GOOntology parseOBO(final Path file) throws IOException
		{
		try(BufferedReader br = IOUtils.openPathForBufferedReading(file)) {
			return parseOBO(br);
			}
		}
	public GOOntology parseOBO(final String uri) throws IOException
		{
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			return parseOBO(br);
			}
		}

	private static class RelationImpl implements GOOntology.Relation
		{
		final TermImpl termFrom;
		final String reltype;
		final TermImpl termTo;
		final private int _hash;
		RelationImpl(final TermImpl termFrom,final String reltype,final TermImpl termTo)
			{
			this.termFrom = termFrom;
			this.reltype = reltype;
			this.termTo = termTo;
			
			final int prime = 31;
			int result = 1;
			result = prime * result + termFrom.hashCode();
			result = prime * result + reltype.hashCode();
			result = prime * result + termTo.hashCode();
			this._hash = result;
			}
		@Override
		public Term getTo() {
			return termTo;
			}
		@Override
		public String getType() {
			return reltype;
			}
		@Override
		public int hashCode() {
			return this._hash;
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			final RelationImpl other = (RelationImpl) obj;
			return (reltype.equals(other.reltype) &&
					termFrom.equals(other.termFrom) &&
					termTo.equals(other.termTo)
					);
			}
		@Override
		public String toString() {
			return this.termFrom.getAcn()+"["+this.termFrom.getName()+"] -- "+
					this.reltype+ " --> "+
					this.termTo.getAcn()+"["+this.termTo.getName()
					;
			}
		}
	
	

	private static class TermImpl implements GOOntology.Term
		{
		final String accession;
		String name;
		String definition = null;
		Set<String> synonyms = null;
		GOOntology.Division division = null;
		//final Set<TermImpl> parents=new HashSet<>();
		//final Set<TermImpl> children=new HashSet<>();
		final Set<RelationImpl> relations = new HashSet<>();
		
		TermImpl(final String accession) {
			this.accession=accession;
			this.name= this.accession;
			this.definition= this.accession;
			}
		
		@Override
		public Set<String> getSynonyms()
			{
			return this.synonyms==null?
					Collections.emptySet():
					Collections.unmodifiableSet(this.synonyms)
					;
			}
		@Override
		public GOOntology.Division getDivision() {
			return this.division;
			}
		
		@Override
		public String getAcn() {
			return accession;
			}
		@Override
		public String getName() {
			return name;
			}
		@Override
		public String getDefinition() {
			return StringUtil.isBlank(definition)?getName():this.definition;
			}

		@Override
		public boolean hasRelations()
			{
			return !this.relations.isEmpty();
			}
		
		@Override
		public Set<GOOntology.Relation> getRelations() {
			return Collections.unmodifiableSet(this.relations);
			}
		
		@Override
		public boolean isDescendantOf(final Term parentNode)
			{
			//LOG.debug(" start"+this+" "+this.relations);
			return _isDescendantOf(parentNode,new HashSet<>());
			}
		
		private boolean _isDescendantOf(final Term parentNode,final Set<Term> seen)
			{
			//LOG.debug("here ?"+this+" "+this.relations);
			if(parentNode==this || getAcn().equals(parentNode.getAcn())) return true;
			seen.add(this);
			for(final RelationImpl rel:this.relations)
				{
				final TermImpl nto = (TermImpl)rel.getTo() ;
				if(seen.contains(nto)) {
					//LOG.debug("twice ?"+this+" "+rel+" "+seen);
					//System.exit(-1);
					continue;
					}
				if(nto._isDescendantOf(parentNode,seen)) return true;
				}
			return false;
			}
		
		@Override
		public int hashCode()
			{
			return this.accession.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			final TermImpl other = (TermImpl) obj;
			return this.accession.equals(other.accession);
			}
		
		@Override
		public int getMinDepth()
			{
			if(!hasRelations()) return 0;
			int d=-1;
			for(final RelationImpl rel:this.relations)
				{
				int d2 = 1 + rel.getTo().getMinDepth();
				if(d==-1 || d>d2)
					{
					d=d2;
					}
				}
			return d;
			}
		
		@Override
		public String toString() {
			return accession;
			}
		}
	
	private static class GOOntologyImpl implements GOOntology {
		private final HashMap<String, TermImpl> acn2term=new HashMap<String, TermImpl>();

		
		GOOntologyImpl()
			{
			}
		
		@Override
		public int size()
			{
			return acn2term.size();
			}
		
		@Override
		public Collection<GOOntology.Term> getTerms() 
			{
			return Collections.unmodifiableCollection(this.acn2term.values());
			}
		
		@Override
		public Iterator<Term> iterator() {
			return getTerms().iterator();
			}
		
		

		
		
		private TermImpl getOrCreate(final String acn) {
			TermImpl t = this.acn2term.get(acn);
			if(t==null) {
				if(!acn.startsWith("GO:")) throw new IllegalArgumentException("Doesn't start with GO:"+acn);
				t = new TermImpl(acn);
				this.acn2term.put(acn,t);
			}
			return t;
			}
		
		@Override
		public Term getTermByAccession(final String s)
			{
			return this.acn2term.get(s);
			}
		
		/** search term by name or accession or synonyms, ignoring case */
		@Override
		public Term getTermByName(final String s)
			{
			final Term t0 = getTermByAccession(s);
			if(t0!=null) return t0;
			for(final TermImpl t:this.acn2term.values())
				{
				if(t.accession.equalsIgnoreCase(s)) return t;
				if(t.name!=null && t.name.equalsIgnoreCase(s)) return t;
				if(t.synonyms!=null)
					{
					for(final String syn:t.synonyms)
						{
						if(syn.equalsIgnoreCase(s)) return t;
						}
					}
				}
			return null;
			}
		}
	
}
