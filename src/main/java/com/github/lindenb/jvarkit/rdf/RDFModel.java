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
package com.github.lindenb.jvarkit.rdf;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.AbstractSet;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Stream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.DC;
import com.github.lindenb.jvarkit.rdf.ns.OWL;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;
import com.github.lindenb.jvarkit.rdf.ns.XSD;

/**
 * Basic Collection of RDF statements
 * @author lindenb
 *
 */
public class RDFModel extends AbstractSet<Statement> implements Set<Statement> {
	private final Set<Statement> statements;
	public RDFModel() {
		statements = new HashSet<>();
		}
	public RDFModel(final Set<Statement> set) {
		statements = new HashSet<>(set);
		}
	@Override
	public boolean addAll(Collection<? extends Statement> c) {
		return statements.addAll(c);
		}
	
	
	public RDFModel addStatement(final Statement e) {
		this.add(e);
		return this;
		}
	
	@Override
	public final boolean add(Statement e) {
		return statements.add(e);
		}

	@Override
	public boolean remove(Object o) {
		return statements.remove(o);
		}
	
	public Set<Resource> getRDFSDescendants(Resource type) {
		return _getRDFSSubClasses(type,RDFS.subClassOf,new HashSet<>());
		}
	
	
	private Set<Resource> _getRDFSSubClasses(Resource type,final Resource predicate,final Set<Resource> set) {
		if(!set.contains(type)) {
			set.add(type);
			findMatching(null, predicate, type).
				map(T->T.getSubject()).
				forEach(S->_getRDFSSubClasses(S,predicate, set));
			}
		return set;
		}
	
	public RDFModel removeStatement(final Resource s, final Resource p, final RDFNode o) {
		remove(new Statement(s,p,o));
		return this;
		}
	
	public RDFModel addStatement(Resource s, Resource p, Date date) {
		return addStatement(new Statement(s,p,new Literal(date)));
		}
	
	public RDFModel addStatement(Resource s, Resource p, RDFNode o) {
		return addStatement(new Statement(s,p,o));
		}
	
	public RDFModel addStatement(Resource s, Resource p, String lit) {
		return addStatement(s,p,new Literal(lit));
		}
	
	public Stream<Statement> findMatching(final Resource s, final Resource p, final RDFNode o) {
		return stream().filter(S->S.match(s, p, o));
		}
	
	public boolean removeMatching(final Resource s, final Resource p, final RDFNode o) {
		return statements.removeIf(S->S.match(s, p, o));
		}
	
	
	@Override
	public Iterator<Statement> iterator() {
		return statements.iterator();
		}
	@Override
	public void clear() {
		statements.clear();
		}

	@Override
	public int size() {
		return statements.size();
		}
	
	
	public void writeXml(final Path outputPath) throws XMLStreamException,IOException {
		try(OutputStream os = IOUtils.openPathForWriting(outputPath)) {
			writeXml(os);
			os.flush();
			}
		}
	
	public void writeXml(final OutputStream out) throws XMLStreamException {
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		final XMLStreamWriter w = xof.createXMLStreamWriter(out, "UTF-8");
		final Map<String,String> ns2prefix = new HashMap<>();
		ns2prefix.put(RDF.NS, RDF.pfx);
		ns2prefix.put(RDFS.NS, RDFS.pfx);
		ns2prefix.put(DC.NS, DC.pfx);
		ns2prefix.put(OWL.NS, OWL.pfx);
		ns2prefix.put(XSD.NS, XSD.pfx);
		
		this.stream().flatMap(S->Arrays.asList(S.getSubject().getNamespaceURI(),S.getPredicate().getNamespaceURI()).stream()).
			forEach(S->{
				if(ns2prefix.containsKey(S)) return;
				final String pfx="_n"+ns2prefix.size();
				ns2prefix.put(S, pfx);
			});
		
		
		w.writeStartDocument("UTF-8", "1.0");
		w.writeStartElement(RDF.pfx,"RDF",RDF.NS);
		for(String ns: ns2prefix.keySet()) {
			w.writeNamespace(ns2prefix.get(ns), ns);
			}
		writeXmlBody(w);
		w.writeEndElement(); //
		w.writeEndDocument();
		}
		
		
	public void writeXmlBody(final XMLStreamWriter w) throws XMLStreamException {
		final Function<Resource, String> toPfx = R->{
			final String pfx;
			try {
				pfx = w.getPrefix(R.getNamespaceURI());
				}
			catch(XMLStreamException err) {
				throw new IllegalStateException(err);
				}
			if(StringUtils.isBlank(pfx)) throw new IllegalStateException("no prefix defined for "+R.getNamespaceURI());
			return pfx;
			};
			
		final Comparator<Statement> compare=(A,B)->A.getSubject().compareTo(A.getSubject());
		
		final EqualIterator<Statement> eq= new EqualIterator<>(
				this.stream().sorted(compare).iterator(),
				compare
				);
		while(eq.hasNext()) {
			final List<Statement> row= eq.next();
			final Statement first = row.get(0);
			final Statement rdfTypeStmt = row.size()==1/* at least 2 stmt */?null:this.findMatching(first.getSubject(), RDF.type, null).filter(S->S.isResource()).findFirst().orElse(null);
			final Resource rdfType = rdfTypeStmt==null?null:rdfTypeStmt.getObject().asResource();
			
			if(rdfType==null) {
				w.writeStartElement(toPfx.apply(RDF.type),"Description", RDF.NS);
				}
			else
				{
				w.writeStartElement( toPfx.apply(rdfType), rdfType.getLocalName(),rdfType.getURI());
				}
			first.getSubject().writeAboutAttribute(w);

			
			
			for(Statement stmt:row) {
				if(stmt.equals(rdfTypeStmt)) continue;
				if(stmt.isLiteral()) {
					final Literal lit = stmt.getObject().asLiteral();
					w.writeStartElement(stmt.getPredicate().getLocalName(), toPfx.apply(stmt.getPredicate()),  stmt.getPredicate().getURI());
					lit.writeRDFXml(w);
					w.writeEndElement();
					}
				else
					{
					w.writeEmptyElement(toPfx.apply(stmt.getPredicate()), stmt.getPredicate().getLocalName(),stmt.getPredicate().getURI());
					stmt.getObject().asResource().writeResourceAttribute(w);
					}			
				}
			
			
			w.writeEndElement();
			}
		}
	}
