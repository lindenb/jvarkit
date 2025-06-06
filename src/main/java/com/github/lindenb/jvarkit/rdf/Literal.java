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

import java.util.Date;
import java.util.Objects;


import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.XSD;

public class Literal implements RDFNode,Comparable<Literal> {
	private final Comparable<?> value;
	public Literal(final String text) {
		this.value = Objects.requireNonNull(text,"text cannot be null");
		}
	public Literal(final long v) {
		this.value = v;
		}
	public Literal(final int v) {
		this.value = v;
		}
	public Literal(final short v) {
		this.value = v;
		}
	public Literal(final float v) {
		this.value = v;
		}
	public Literal(final double v) {
		this.value = v;
		}
	public Literal(final boolean v) {
		this.value = v;
		}
	public Literal(final Date v) {
		this.value = v;
		}
	
	@Override
	public final boolean isResource() {
		return false;
		}
	@Override
	public final boolean isLiteral() {
		return true;
		}
	
	/** return the internal java object */
	public Comparable<?> getValue() {
		return value;
		}
	
	public boolean isString() {
		return value instanceof String;
		}
	public boolean isDate() {
		return value instanceof Date;
		}
	
	public String getDatatypeURI() {
		if(isString()) return XSD.NS+"string";
		if(isDate()) return XSD.NS+"date";
		throw new IllegalStateException("not a string");
		}
	
	/** return the internal object as string whatever is its class */
	public String getLiteralAsString() {
		return String.valueOf(getValue());
		}
	
	
	public String getString() {
		if(isString()) {
			return String.class.cast(value);
			}
		throw new IllegalStateException("not a string. Use getLiteralAsString ?");
		}
	
	public Date getDate() {
		if(isDate()) {
			return Date.class.cast(getValue());
			}
		throw new IllegalStateException("not a Date. Use getLiteralAsString ?");
		}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public int compareTo(final Literal o) {
		if(o==this) return 0;
		final Class<?> C1 = this.getValue().getClass();
		final Class<?> C2 = o.getValue().getClass();
		if(!C1.equals(C2)) {
			return C1.getName().compareTo(C2.getName());
			}
		final Comparable v1 = this.getValue();
		final Comparable v2 = o.getValue();
		return v1.compareTo(v2);
		}
	
	@Override
	public int hashCode() {
		return Objects.hash(value);
		}
	
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null || !(obj instanceof Literal))
			return false;
		final Literal other = (Literal) obj;
		return Objects.equals(getValue(), other.getValue());
		}
	
	@Override
	public String toString() {
		return new StringBuilder("\"").
				append(StringUtils.escapeC(getLiteralAsString())).
				append("\"").
				toString();
		}
	}
