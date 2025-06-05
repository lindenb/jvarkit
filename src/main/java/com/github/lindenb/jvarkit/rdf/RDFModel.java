package com.github.lindenb.jvarkit.rdf;

import java.io.OutputStream;
import java.util.AbstractSet;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
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

import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.rdf.ns.DC;
import com.github.lindenb.jvarkit.rdf.ns.OWL;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;
import com.github.lindenb.jvarkit.rdf.ns.XSD;

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
	
	
	public RDFModel addStatement(Statement e) {
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
	
	public RDFModel removeStatement(Resource s, Resource p, RDFNode o) {
		remove(new Statement(s,p,o));
		return this;
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
	
	public void writeXml(OutputStream out) throws XMLStreamException {
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		final XMLStreamWriter w = xof.createXMLStreamWriter(out, "UTF-8");
		final Map<String,String> ns2prefix = new HashMap<>();
		final Function<Resource, String> toPfx = R->ns2prefix.get(R.getNamespaceURI());
		ns2prefix.put(RDF.NS, RDF.pfx);
		ns2prefix.put(RDFS.NS, RDFS.pfx);
		ns2prefix.put(DC.NS, DC.pfx);
		ns2prefix.put(OWL.NS, OWL.pfx);
		ns2prefix.put(XSD.NS, XSD.pfx);
		
		this.stream().flatMap(S->Arrays.asList(S.getSubject().getNamespaceURI(),S.getPredicate().getNamespaceURI()).stream()).
			forEach(S->{
				if(ns2prefix.containsKey(S)) return;
				String pfx="_n"+ns2prefix.size();
				ns2prefix.put(S, pfx);
			});
		
		
		w.writeStartDocument("UTF-8", "1.0");
		w.writeStartElement("rdf:RDF");
		
		for(String ns: ns2prefix.keySet()) {
			w.writeNamespace(ns2prefix.get(ns), ns);
			}
		
		
		final Comparator<Statement> compare=(A,B)->A.getSubject().compareTo(A.getSubject());
		
		EqualIterator<Statement> eq= new EqualIterator<>(
				this.stream().sorted(compare).iterator(),
				compare
				);
		while(eq.hasNext()) {
			final List<Statement> row= eq.next();
			final Statement first = row.get(0);
			final Statement rdfTypeStmt = row.size()==1/* at least 2 stmt */?null:this.findMatching(first.getSubject(), Resource.RDF_type, null).filter(S->S.isResource()).findFirst().orElse(null);
			final Resource rdfType = rdfTypeStmt==null?null:rdfTypeStmt.getObject().asResource();
			
			if(rdfType==null) {
				w.writeStartElement("rdf","Description", RDF.NS);
				}
			else
				{
				w.writeStartElement( toPfx.apply(rdfType), rdfType.getLocalName(),rdfType.getURI());
				}
			w.writeAttribute("rdf", RDF.NS, "about", first.getSubject().getURI());

			
			
			for(Statement stmt:row) {
				if(stmt.equals(rdfTypeStmt)) continue;
				if(stmt.isLiteral()) {
					final Literal lit = stmt.getObject().asLiteral();
					w.writeStartElement(stmt.getPredicate().getLocalName(), toPfx.apply(stmt.getPredicate()),  stmt.getPredicate().getURI());
					if(!lit.isString()) {
						w.writeAttribute("rdf", RDF.NS, "dataType",lit.getDatatypeURI());
						}
					w.writeCharacters(lit.getString());
					w.writeEndElement();
					}
				else
					{
					w.writeEmptyElement(toPfx.apply(stmt.getPredicate()), stmt.getPredicate().getLocalName(),stmt.getPredicate().getURI());
					w.writeAttribute("rdf", RDF.NS, "resource", stmt.getObject().asResource().getURI());
					}
			
				}
			
			
			w.writeEndElement();
			}
		
		
		w.writeEndElement(); //
		w.writeEndDocument();
		}
	}
