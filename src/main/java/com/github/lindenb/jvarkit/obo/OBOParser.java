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
package com.github.lindenb.jvarkit.obo;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;

/** direct graph for Open Biological and Biomedical Ontology  (OBO) */
public class OBOParser {
	private static final Logger LOG = Logger.build(OBOParser.class).make();

		private static class OBOntologyImpl implements OBOntology {
			private final Map<String,OBOntologyImpl.TermImpl> acn2term = new HashMap<>();
			
			private class TermImpl implements OBOntology.Term 	{
				private final String id;
				private final String label;
				private final int _hash;
				final Set<TermImpl> parents=new HashSet<>();
				final Set<TermImpl> children=new HashSet<>();
		
				TermImpl(final String id,final String label) {
					this.id = id;
					this.label = label;
					this._hash = id.hashCode();
					}
				@Override
				public String getLabel() {
					return this.label;
					}
				@Override
				public String getAcn() {
					return this.id;
					}
				@Override
				public int hashCode() {
					return _hash;
					}
				@Override
				public boolean equals(final Object obj) {
					if(obj==this) return true;
					if(obj==null || !(obj instanceof TermImpl)) return false;
					return this._hash==obj.hashCode() && getAcn().equals(TermImpl.class.cast(obj).getAcn());
					}
				@Override
				public String toString() {
					return getAcn();
					}
				@Override
				public Set<Term> getParents() {
					return Collections.unmodifiableSet(this.parents);
					}
				@Override
				public Set<Term> getChildren() {
					return Collections.unmodifiableSet(this.children);
					}
				
				private Set<Term> visitParents(final Set<Term> set) {
					if(set.contains(this)) return set;
					set.add(this);
					for(final TermImpl t:this.parents) t.visitParents(set);
					return set;
					}
				private Set<Term> visitChildren(final Set<Term> set) {
					if(set.contains(this)) return set;
					set.add(this);
					for(final TermImpl t:this.children) t.visitParents(set);
					return set;
					}
				
				@Override
				public Set<Term> getAllAncestors() {
					return visitParents(new HashSet<>());
					}
				@Override
				public Set<Term> getAllDescendant() {
					return visitChildren(new HashSet<>());
					}
				
				private boolean searchAncestor(final Set<Term> set,final Term node) {
					if(set.contains(node)) return false;
					if(node.equals(this)) return true;
					set.add(this);
					for(final TermImpl t:this.parents) if(t.searchAncestor(set,node)) return true;
					return false;
					}
				
				@Override
				public boolean isAncestorOf(final Term node) {
					return searchAncestor(new HashSet<>(),node);
					}
				
				@Override
				public final boolean isDescendantOf(final Term node) {
					return node.isAncestorOf(this);
					}
				}
			
			@Override
			public Stream<Term> stream() {
				return this.acn2term.
						values().
						stream().
						map(T->Term.class.cast(T));
				}
			@Override
			public Term findTermByAcn(final String id) {
				return this.acn2term.get(id);
				}
			@Override
			public Iterator<Term> iterator() {
				return stream().iterator();
				}
			}

		private final String OWL="http://www.w3.org/2002/07/owl#";
		private final String RDF="http://www.w3.org/1999/02/22-rdf-syntax-ns#";
		private final String RDFS="http://www.w3.org/2002/07/owl#";
		private final String  oboInOwl="http://www.geneontology.org/formats/oboInOwl#";
		private final QName rdf_about=  new QName(RDF, "about");
		private final QName rdf_resource=  new QName(RDF, "about");
		private Map<String,OBOntologyImpl.TermImpl> uri2term= new HashMap<>();
		private List<Object> isAList= new ArrayList<>();
		
		private void skip(final XMLEventReader r) throws XMLStreamException
			{
			while(r.hasNext())
				{
				final XMLEvent evt = r.nextEvent();
				if(evt.isEndElement()) break;
				else if(evt.isStartElement())
					{
					skip(r);
					}
				}
			}
	

		
		private boolean isA(final QName qName,final String ns,final String localName) {
			return qName.getNamespaceURI().equals(ns) && qName.getLocalPart().equals(localName);
			}
		
		private void parseOWLClass(
				final OBOntologyImpl ontology,
				final StartElement SE,
				final XMLEventReader r
				) throws IOException,XMLStreamException {
			
			Attribute att = SE.getAttributeByName(rdf_about);
			String label=null;
			String id=null;
			String about = (att==null?null:att.getValue());
			boolean deprecated=false;
			final Set<String> subClassesOf=new HashSet<>(); 
			while(r.hasNext())
				{
				final XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					final StartElement E=evt.asStartElement();
					final QName qN=E.getName();
					if(isA(qN,RDFS,"label")) {
						label = r.getElementText(); 
						}
					else if(isA(qN,OWL,"deprecated")) {
						deprecated = r.getElementText().equals("true"); 
						}
					else if(isA(qN,oboInOwl,"id")) {
						id = r.getElementText();
						}
					else if(isA(qN,RDFS,"subClassOf")) {
						att = SE.getAttributeByName(rdf_resource);
						if(att!=null)
							{
							subClassesOf.add(att.getValue());
							}
						}
					else
						{
						skip(r);
						}
					}
				else if(evt.isEndElement()) {
					break;
					}
				}
			if(StringUtils.isBlank(id)) return;
			if(StringUtils.isBlank(label)) return;
			if(StringUtils.isBlank(about)) return;
			if(deprecated) return;
			if(this.uri2term.containsKey(about)) {
				throw new XMLStreamException("duplicate owl:Class "+about,SE.getLocation());
				}
			final OBOntologyImpl.TermImpl term =ontology.new TermImpl(id,label);
			this.uri2term.put(about, term);
			for(final String x:subClassesOf) {
				this.isAList.add(term);
				this.isAList.add(x);
				}
			}
		
		public OBOntology parseOwl(final Reader xmlIn) throws IOException,XMLStreamException
			{
			XMLEventReader r=null;
			try
				{
				final OBOntologyImpl ontology=new OBOntologyImpl();
				final XMLInputFactory fact=XMLInputFactory.newFactory();
				fact.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			    r=fact.createXMLEventReader(xmlIn);
				while(r.hasNext())
					{
					final XMLEvent evt=r.nextEvent();
					if(evt.isStartElement())
						{
						final StartElement E=evt.asStartElement();
						final QName qN=E.getName();
						if(qN.getLocalPart().equals("Class") && 
							OWL.equals(qN.getNamespaceURI()))
							{
							this.parseOWLClass(ontology,E,r);
							}
						}
					}
				r.close();
				for(final OBOntologyImpl.TermImpl t:this.uri2term.values()) {
					ontology.acn2term.put(t.getAcn(), t);
					}
				
				for(int i=0;i+1< this.isAList.size();i+=2)
					{
					final OBOntologyImpl.TermImpl term = OBOntologyImpl.TermImpl.class.cast(this.isAList.get(i+0));
					final String parentId = String.class.cast(this.isAList.get(i+1));
					final OBOntologyImpl.TermImpl parent = this.uri2term.get(parentId);
					if(parent==null) {
						System.err.println("Cannot find uri :"+parentId);
						continue;
						}
					if(term.equals(parent))
						{
						LOG.warning("self/self relation ??"+parentId);
						continue;
						}
					parent.children.add(term);
					term.parents.add(parent);
					}				
				r.close();
				return ontology;
				}
			finally
				{
				this.uri2term.clear();
				this.isAList.clear();
				CloserUtil.close(r);
				}
			}

		public OBOntology parseOwl(final Path path) throws IOException,XMLStreamException
			{
			try(Reader r=IOUtils.openFileForReader(path)) {
				return parseOwl(r);
				}
			}
		
		/** parse ontology using file extension. '.obo' or 'obo.gz' are interpretted as OBO. Default is RDF/XML */
		public OBOntology parseByExtension(final Path path) throws IOException,XMLStreamException {
			String suff  = path.getFileName().toString();
			if(suff.endsWith(".obo") || suff.endsWith(".obo.gz")) {
				return parseObo(path);
				}
			else
				{
				return parseOwl(path);
				}
			}
		
		
		public OBOntology parseObo(final Path path) throws IOException,XMLStreamException
			{
			try(BufferedReader r=IOUtils.openPathForBufferedReading(path)) {
				return parseObo(r);
				}
			}

		
		public OBOntology parseObo(final BufferedReader br) throws IOException,XMLStreamException
			{
			try
				{
				final OBOntologyImpl ontology=new OBOntologyImpl();
				final Iterator<String> iter1 = IOUtils.toLineIterator(br);
				final OBOIterator iter = new OBOIterator(iter1);
				while(iter.hasNext())
					{
					final Map<String,Set<String>> hash = iter.next();
					
					Set<String> set = hash.get("name");
					if(set.size()!=1) continue;
					final String label= set.iterator().next();
					
					set = hash.get("id");
					if(set.size()!=1) continue;
					final String id= set.iterator().next();
					
					final OBOntologyImpl.TermImpl term =ontology.new TermImpl(id,label);
					ontology.acn2term.put(id, term);
					
					set = hash.get("is_a");
					
					
					for(final String x:set) {
						this.isAList.add(term);
						this.isAList.add(x);
						}

					}
				
				br.close();
				
				
				for(int i=0;i+1< this.isAList.size();i+=2)
					{
					final OBOntologyImpl.TermImpl term = OBOntologyImpl.TermImpl.class.cast(this.isAList.get(i+0));
					final String parentId = String.class.cast(this.isAList.get(i+1));
					final OBOntologyImpl.TermImpl parent = ontology.acn2term.get(parentId);
					if(parent==null) {
						LOG.warning("Cannot find uri :"+parentId);
						continue;
						}
					if(term.equals(parent))
						{
						LOG.warning("self/self relation ??"+parentId);
						continue;
						}
					parent.children.add(term);
					term.parents.add(parent);
					}				
				return ontology;
				}
			finally
				{
				this.uri2term.clear();
				this.isAList.clear();
				}
			}

		
	
		private static class OBOIterator extends AbstractIterator<Map<String,Set<String>>>
			{
			final Map<String,Set<String>> buffer= new HashMap<>();
			final Iterator<String> delegate;
			boolean in_term=false;
			
			OBOIterator(final Iterator<String> delegate) {
				this.delegate=delegate;
				}
			
			private Map<String,Set<String>> dump() {
				final Map<String,Set<String>> L = new HashMap<>(this.buffer);
				buffer.clear();
				in_term=false;
				return L;
				}
			@Override
			protected Map<String,Set<String>> advance() {
				for(;;)
					{
					final String line = this.delegate.hasNext()?
							this.delegate.next():
							null;
					
					if(StringUtils.isBlank(line))
						{
						if(in_term && !this.buffer.isEmpty()) return dump();
						if(line==null) return null;
						this.buffer.clear();
						this.in_term=false;
						}
					else if(line.trim().equals("[Term]"))
						{
						if(in_term && !this.buffer.isEmpty()) return dump();
						this.buffer.clear();
						this.in_term=true;
						}
					else if(this.in_term)
						{
						int idx= line.indexOf(":");
						if(idx==-1) continue;
						String key= StringUtils.substringBefore(line, ":").trim();
						String value = StringUtils.substringAfter(line, ":").trim();
						Set<String> set = this.buffer.get(key);
						if(set==null) {
							set = new HashSet<>();
							this.buffer.put(key,set);
							}
						if(key.equals("is_a") && value.contains("!"))
							{
							value= StringUtils.substringBefore(value,"!").trim();
							}
						set.add(value);
						}
					}
				}
			}
	
}
