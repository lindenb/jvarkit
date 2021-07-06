/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;

public class UKBiobankOntology {
	private static final Logger LOG = Logger.build(UKBiobankOntology.class).make();

	/** a term in the ontology */
	public interface Term {
		public String getCoding();
		public String getMeaning();
		public Term getParent();
		public boolean isSelectable();
		public Set<Term> getAllChildrenIncludingSelf();
		//public Set<Term> getAllParentsIncludingSelf();
		public int getDepth();
	}
	
	private class TermImpl implements Term {
		String coding="";
		String meaning="";
		String node_id="";
		String parent_id="";
		boolean selectable=true;
		private Set<Term> _children_cache = null;
		
		@Override
		public int hashCode() {
			return node_id.hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			if(obj==null || !(obj instanceof TermImpl)) return false;
			final TermImpl other = TermImpl.class.cast(obj);
			return coding.equals(other.coding) &&
					meaning.equals(other.meaning) &&
					node_id.equals(other.node_id) &&
					parent_id.equals(other.parent_id) &&
					selectable == other.selectable
					;
			}
		@Override
		public boolean isSelectable() {
			return this.selectable;
			}
		@Override
		public String getCoding() {
			return coding;
			}
		@Override
		public String getMeaning() {
			return meaning;
			}
		@Override
		public Term getParent() {
			return UKBiobankOntology.this.nodeid2term.get(this.parent_id);
			}
		@Override
		public Set<Term> getAllChildrenIncludingSelf() {
			if(this._children_cache == null) {
				this._children_cache = new HashSet<>();
				scanForChildren(this.node_id,this._children_cache);
				}
			return this._children_cache;
			}
		/*
		public Set<Term> getAllParentsIncludingSelf() {
			final Set<Term> set = new HashSet<>();
			scanForParents(this.node_id,set);
			return set;
			}*/
		@Override
		public int getDepth() {
			Term parent = getParent();
			return (parent==null?0:1+parent.getDepth());
		}
		
		@Override
		public String toString() {
			return getCoding()+"{"+getMeaning()+"}";
			}
		}
	
	private final Map<String, TermImpl> nodeid2term = new HashMap<>();
	private final Map<String, TermImpl> coding2term = new HashMap<>();
	private Path sourcePath=null;

	private UKBiobankOntology() {
	}
	
	private void scanForChildren(final String nodeId,final Set<Term> found) {
			if(!this.nodeid2term.containsKey(nodeId)) {
				LOG.warn("cannot find node_id="+nodeId+" in "+ UKBiobankOntology.this.sourcePath);
				return;
			}
			final TermImpl parent = this.nodeid2term.get(nodeId);
			found.add(parent);
			for(TermImpl child:this.nodeid2term.values()) {
				if(child.parent_id.equals(nodeId)) {
					scanForChildren(child.node_id,found);
				}
			}
		}
	/*
	private void scanForParents(String nodeId,final Set<Term> found) {
		if(!this.nodeid2term.containsKey(nodeId)) {
			LOG.warn("cannot find node_id="+nodeId+" in "+ UKBiobankOntology.this.sourcePath);
			return;
			}
		Term node = this.nodeid2term.get(nodeId);
		found.add(node);
		for(;;) {
			node = node.getParent();
			if(node==null) break;
			found.add(node);
		}
	}*/

	public Term findTermByCoding(final String coding) {
		return this.coding2term.get(coding);
		}
	
	
	public Path getSource() {
		return sourcePath;
		}
	
	@Override
	public String toString() {
		return getSource().toString();
		}
	
	public static UKBiobankOntology load(final Path path) throws IOException {
		final String[] expect_cols=new String[] {"coding","meaning","node_id","parent_id","selectable"};
		final UKBiobankOntology ontology = new UKBiobankOntology();
		ontology.sourcePath = path;
		try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
			String line = br.readLine();
			if(line==null) throw new IOException("cannot read first line of "+path);
			String[] tokens = CharSplitter.TAB.split(line);
			if(tokens.length!=expect_cols.length) throw new JvarkitException.TokenErrors(expect_cols.length, tokens);
			for(int x=0;x< expect_cols.length;++x) {
				if(!expect_cols[x].equals(tokens[x])) throw new IOException("expected "+expect_cols[x]+" but got "+tokens[x]+" in "+line);
			}
			
			while((line=br.readLine())!=null) {
				tokens = CharSplitter.TAB.split(line);
				if(tokens.length!=expect_cols.length) throw new JvarkitException.TokenErrors(expect_cols.length, tokens);
				final TermImpl term = ontology.new TermImpl();
				term.coding = tokens[0];
				term.meaning = tokens[1];
				term.node_id = tokens[2];
				term.parent_id = tokens[3];
				term.selectable = tokens[4].equals("Y");
				if(ontology.nodeid2term.containsKey(term.node_id)) {
					throw new IOException("duplicate node id for "+line+" in "+path);
					}
				if(ontology.coding2term.containsKey(term.coding)) {
					throw new IOException("duplicate coding for "+line+" in "+path);
					}
				ontology.nodeid2term.put(term.node_id, term);
				ontology.coding2term.put(term.coding, term);
			}
		}
		return ontology;
	}
	
}
