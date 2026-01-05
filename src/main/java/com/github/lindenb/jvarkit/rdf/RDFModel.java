/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import java.io.Writer;
import java.nio.file.Path;
import java.util.AbstractSet;
import java.util.ArrayList;
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
import com.github.lindenb.jvarkit.util.Algorithm;

/**
 * Basic Collection of RDF statements
 * @author lindenb
 *
 */
public class RDFModel extends AbstractSet<Statement> implements Set<Statement> {
	private final List<Statement> sorted_statements = new ArrayList<>();
	
	
	private final static Algorithm<Statement, Resource> SUBJECT_SORTER = new Algorithm<>(
			STMT->STMT.getSubject(),
			(A,B)-> A.compareTo(B)
			);
	
	
	private final static Algorithm<Statement, Statement> STMT_SORTER = new Algorithm<>(
			STMT->STMT,
			(A,B)->{
				int i = SUBJECT_SORTER.getComparator().compare(A.getSubject(),B.getSubject());
				if(i!=0) return i;
				i = A.getPredicate().compareTo(B.getPredicate());
				if(i!=0) return i;
				if(A.getObject().isResource()) {
					if(B.isResource()) {
						return A.getObject().asResource().compareTo(B.getObject().asResource());
						}
					else if(B.isLiteral())  {
						return -1;
						}
					else
						{
						throw new IllegalArgumentException();
						}
					}
				else if(A.getObject().isLiteral()) {
					if(B.isResource()) {
						return 1;
						}
					else if(B.isLiteral())  {
						return A.getObject().asLiteral().compareTo(B.getObject().asLiteral());
						}
					else
						{
						throw new IllegalArgumentException();
						}
					}
				else
					{
					throw new IllegalArgumentException();
					}
				}
			);
	
	public RDFModel() {
		}
	public RDFModel(final Set<Statement> set) {
		this.addAll(set);
		}
	
	@Override
	protected RDFModel clone() {
		final RDFModel m=new RDFModel();
		m.sorted_statements.addAll(this.sorted_statements);
		return m;
		}
	
	@Override
	public boolean addAll(final Collection<? extends Statement> c) {
		for(Statement stmt: c) {
			add(stmt);
			}
		return true;
		}
	
	
	public RDFModel addStatement(final Statement e) {
		this.add(e);
		return this;
		}
	
	@Override
	public final boolean add(final Statement stmt) {
		final int idx = STMT_SORTER.lower_bound(sorted_statements,stmt);
		if(idx < sorted_statements.size()) {
			if(sorted_statements.get(idx).equals(stmt)) return false;
			}
		sorted_statements.add(idx,stmt);
		return true;
		}

	@Override
	public boolean contains(final Object o) {
		if(o==null || !(o instanceof Statement)) return false;
		final Statement stmt = Statement.class.cast(o);
		return containsStatement(stmt);
		}
	
	public boolean containsStatement(Resource subject, Resource predicate, RDFNode obj) {
		return containsStatement(new Statement(subject,predicate,obj));
		}

	
	public boolean containsStatement(final Statement stmt) {
		final int idx = STMT_SORTER.lower_bound(sorted_statements,stmt);
		if(idx < sorted_statements.size()) {
			if(sorted_statements.get(idx).equals(stmt)) {
				return true;
				}
			}
		return false;
		}
	
	
	@Override
	public boolean remove(final Object o) {
		if(o==null || !(o instanceof Statement)) return false;
		final Statement stmt = Statement.class.cast(o);
		final int idx = STMT_SORTER.lower_bound(sorted_statements,stmt);
		if(idx < sorted_statements.size()) {
			if(sorted_statements.get(idx).equals(stmt)) {
				sorted_statements.remove(idx);
				return true;
				}
			}
		return false;
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
	
	public Stream<Statement> findMatching(final Resource subject, final Resource p, final RDFNode o) {
		if(subject!=null) {
			final int[] idxs = SUBJECT_SORTER.equal_range(this.sorted_statements, subject);
			return this.sorted_statements.subList(idxs[0],idxs[1]).stream().filter(S->S.match(subject, p, o));
			}
		return stream().filter(S->S.match(subject, p, o));
		}
	
	public boolean removeMatching(final Resource subject, final Resource p, final RDFNode o) {
		if(subject!=null) {
			final int[] idxs = SUBJECT_SORTER.equal_range(this.sorted_statements, subject);
			return this.sorted_statements.subList(idxs[0],idxs[1]).removeIf(S->S.match(subject, p, o));
			}
		return this.sorted_statements.removeIf(S->S.match(subject, p, o));
		}
	
	
	@Override
	public Iterator<Statement> iterator() {
		return sorted_statements.iterator();
		}
	@Override
	public void clear() {
		sorted_statements.clear();
		}

	@Override
	public int size() {
		return sorted_statements.size();
		}
	
	/** export model to XML+RDF */
	public void writeXml(final Path outputPath) throws XMLStreamException,IOException {
		try(OutputStream os = IOUtils.openPathForWriting(outputPath)) {
			writeXml(os);
			os.flush();
			}
		}
	
	/** create Mapping NS->PREFIX */
	private Map<String,String> getNamespaceToPrefixMap() {
		final Map<String,String> ns2prefix = new HashMap<>();
		ns2prefix.put(RDF.NS, RDF.pfx);
		ns2prefix.put(RDFS.NS, RDFS.pfx);
		ns2prefix.put(DC.NS, DC.pfx);
		ns2prefix.put(OWL.NS, OWL.pfx);
		ns2prefix.put(XSD.NS, XSD.pfx);
		ns2prefix.put("http://www.w3.org/2004/02/skos/core#", "skos");
		ns2prefix.put("http://purl.org/dc/terms/", "dcterms");
		ns2prefix.put("http://xmlns.com/foaf/0.1/", "foaf");

				
		this.stream().flatMap(S->Arrays.asList(S.getSubject().getNamespaceURI(),S.getPredicate().getNamespaceURI()).stream()).
			forEach(S->{
				if(ns2prefix.containsKey(S)) return;
				final String pfx="_n"+ns2prefix.size();
				ns2prefix.put(S, pfx);
			});
		return ns2prefix;
		}
	
	/** export RDF as Turtle format */
	public void writeTurtle(final Path outputPath) throws IOException {
		try(Writer os = IOUtils.openPathForPrintWriter(outputPath)) {
			writeTurtle(os);
			os.flush();
			}
		}
	
	/** export RDF as Turtle format */
	public void writeTurtle(final Writer out) throws IOException {
		final Map<String,String> ns2prefix = getNamespaceToPrefixMap();
		for(String ns: ns2prefix.keySet()) {
			out.append("@prefix ").
			append(ns2prefix.get(ns)).
				append(": <").append(ns).
				append("> .\n");
			}
		out.append("\n");
		
		final Function<Resource,String> rsrc2uri=RSRC->{
			final String pfx = ns2prefix.get(RSRC.getNamespaceURI());
			if(!StringUtils.isBlank(pfx)) return pfx+":"+RSRC.getLocalName();
			return "<"+RSRC.getURI()+">";
			};
		final Comparator<Statement> compare=(A,B)->A.getSubject().compareTo(A.getSubject());
		
		final EqualIterator<Statement> eq= new EqualIterator<>(iterator(), compare );
		while(eq.hasNext()) {
			List<Statement> row = eq.next();
			final Statement first = row.get(0);
			out.write(rsrc2uri.apply(first.getSubject()));
			out.write("\n");
			for(int i=0;i< row.size();i++) {
				out.write("    ");
				final Statement stmt = row.get(i);
				if(stmt.getPredicate().equals(RDF.type)) {
					out.write("a");
					}
				else {
					out.write(rsrc2uri.apply(stmt.getPredicate()));
					}
			
				out.write(" ");
				if(stmt.isLiteral()) {
					final Literal lit = stmt.getObject().asLiteral();
					out.write(StringUtils.doubleQuote(lit.getLiteralAsString()));
					if(lit.hasLang()) {
						out.write("@");
						out.write(StringUtils.doubleQuote(lit.getLang()));
						}
				} else {
					out.write(rsrc2uri.apply(stmt.getObject().asResource()));
					}
				out.write(i+1==row.size()?" .\n":" ;\n");
				}
			out.write("\n");
			}
				
		out.flush();
		}
	
	public void writeXml(final OutputStream out) throws XMLStreamException {
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		final XMLStreamWriter w = xof.createXMLStreamWriter(out, "UTF-8");
		final Map<String,String> ns2prefix = getNamespaceToPrefixMap();
		
		
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
		
		final EqualIterator<Statement> eq= new EqualIterator<>(iterator(), compare );
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
			if(!first.getSubject().isAnon()) {
				first.getSubject().writeAboutAttribute(w);
				}
			
			
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
