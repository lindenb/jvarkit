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

import java.util.Date;
import java.util.Objects;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;
/**
 * RDF Statement(S,P,O)
 *
 */
public class Statement {
	private final Resource subject;
	private final Resource property;
	private final RDFNode object;
	
	public Statement(Resource subject,Resource property,RDFNode object) {
		this.subject = Objects.requireNonNull(subject);
		this.property = Objects.requireNonNull(property);
		this.object = Objects.requireNonNull(object);
		}
	public Statement(Resource subject,Resource property,final Date object) {
		this(subject, property, new Literal(object));
		}
	public Statement(Resource subject,Resource property,final String object) {
		this(subject, property, new Literal(object));
		}
	public Statement(Resource subject,Resource property,final int object) {
		this(subject, property, new Literal(object));
		}
	public Statement(Resource subject,Resource property,final long object) {
		this(subject, property, new Literal(object));
		}
	public Statement(Resource subject,Resource property,final float object) {
		this(subject, property, new Literal(object));
		}
	public Statement(Resource subject,Resource property,final double object) {
		this(subject, property, new Literal(object));
		}
	
	public boolean match(final Resource subject,final Resource property,final RDFNode object) {
		return (subject==null || subject.equals(this.subject)) &&
				(property==null || property.equals(this.property)) &&
				(object==null || object.equals(this.object))
				;
		}
	
	
	
	public boolean isLiteral() {
		return getObject().isLiteral();
		}
	public boolean isResource() {
		return getObject().isResource();
		}
	
	public Resource getSubject() {
		return subject;
		}
	public Resource getPredicate() {
		return property;
		}
	public RDFNode getObject() {
		return object;
		}
	@Override
	public int hashCode() {
		return Objects.hash(subject,property,object);
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Statement)) return false;
		final Statement stmt=Statement.class.cast(obj);
		return this.getSubject().equals(stmt.getSubject()) &&
				 this.getPredicate().equals(stmt.getPredicate()) &&
				 this.getObject().equals(stmt.getObject())
				 ;
		}
	
	/** write this as &lt;rdf:Statement/&gt; */
	public void writeRDFStatement(XMLStreamWriter w) throws XMLStreamException {
		final String pfx= StringUtils.ifBlank(w.getPrefix(RDFS.NS),RDF.pfx);
		w.writeStartElement(pfx, "Statement", RDF.NS);
		
		w.writeEmptyElement(pfx, "subject", RDF.NS);
		this.getSubject().writeResourceAttribute(w);
		
		w.writeEmptyElement(pfx, "predicate", RDF.NS);
		this.getPredicate().writeResourceAttribute(w);
		
		
		if(isLiteral()) {
			final Literal lit = getObject().asLiteral();
			w.writeStartElement(pfx, "object", RDF.NS);
			lit.writeRDFXml(w);
			w.writeEndElement();
			}
		else
			{
			w.writeEmptyElement(pfx, "object", RDF.NS);
			this.getObject().asResource().writeResourceAttribute(w);
			}
		w.writeEndElement();
		}
	
	@Override
	public String toString() {
		return new StringBuilder(this.subject.toString()).append(" ").
				append(this.property.toString()).append(" ").
				append(this.object.toString()).append(" .").
				toString();
		}
	}
