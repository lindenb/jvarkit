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

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Date;
import java.util.Objects;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.XSD;

/**
 * 
 * RDF literal
 *
 */
public class Literal implements RDFNode {
	private final Comparable<?> value;
	private final String lang ;
	public Literal(final String text,final String lang) {
		this.value = Objects.requireNonNull(text,"text cannot be null");
		this.lang  =lang;
		}
	public Literal(final String text) {
		this(text,null);
		}
	public Literal(final long v) {
		this.value = v;
		this.lang  = null;
		}
	public Literal(final int v) {
		this.value = v;
		this.lang  = null;
		}
	
	public Literal(final short v) {
		this.value = v;
		this.lang  = null;
		}
	public Literal(final byte v) {
		this.value = v;
		this.lang  = null;
		}
		
	public Literal(final float v) {
		this.value = v;
		this.lang  = null;
		}
	public Literal(final double v) {
		this.value = v;
		this.lang  = null;
		}
	public Literal(final boolean v) {
		this.value = v;
		this.lang  = null;
		}
	
	public Literal(final BigInteger v) {
		this.value = Objects.requireNonNull(v);
		this.lang  = null;
		}
	
	public Literal(final BigDecimal v) {
		this.value = Objects.requireNonNull(v);
		this.lang  = null;
		}
	
	public Literal(final Date v) {
		this.value = v;
		this.lang  = null;
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
	
	/** return wether this literal is associated to a lang */
	public boolean hasLang() {
		return isString() && !StringUtils.isBlank(this.lang);
	}
	
	/** return lang, may be empty  */
	public String getLang() {
		if(!isString()) throw new IllegalStateException("no a string");
		return StringUtils.ifBlank(this.lang,"");
		}
	
	public boolean isString() {
		return value instanceof String;
		}
	public boolean isDate() {
		return value instanceof Date;
		}
	
	public boolean isNumber() {
		return value instanceof Number;
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
	
	public long getLong() {
		if(isNumber()) return getNumber().longValue();
		throw new IllegalStateException("not along");
		}
	
	public Number getNumber() {
		if(isNumber()) {
			return Number.class.cast(getValue());
			}
		throw new IllegalStateException("not a Number. Use getLiteralAsString ?");
		}
	
	public boolean isInteger() {
		//NON if(this.value instanceof Boolean) return true;
		if(this.value instanceof Byte) return true;
		if(this.value instanceof Short) return true;
		if(this.value instanceof Integer) return true;
		if(this.value instanceof Long) return true;
		if(this.value instanceof BigInteger) return true;
		return false;
		}
	
	public BigInteger asBigInteger() {
		//NON if(this.value instanceof Boolean) return BigInteger.valueOf(Boolean.class.cast(this.value).equals(Boolean.TRUE)?1:0);
		if(this.value instanceof Byte) return BigInteger.valueOf(Byte.class.cast(this.value));
		if(this.value instanceof Short) return BigInteger.valueOf(Short.class.cast(this.value));
		if(this.value instanceof Integer) return BigInteger.valueOf(Integer.class.cast(this.value));
		if(this.value instanceof Long) return BigInteger.valueOf(Long.class.cast(this.value));
		if(this.value instanceof BigInteger) return BigInteger.class.cast(this.value);
		throw new IllegalStateException("not an Integer");
		}
	
	public boolean isFloating() {
		if(isInteger()) return true;
		if(this.value instanceof Float) return true;
		if(this.value instanceof Double) return true;
		if(this.value instanceof BigDecimal) return false;
		return false;
		}
	
	
	public BigDecimal asBigDecimal() {
		if(isInteger()) {
			return new BigDecimal(asBigInteger());
			}
		if(this.value instanceof Float) return BigDecimal.valueOf(Float.class.cast(this.value));
		if(this.value instanceof Double) return BigDecimal.valueOf(Double.class.cast(this.value));
		if(this.value instanceof BigDecimal) return BigDecimal.class.cast(this.value);
		throw new IllegalStateException("not an BigDecimal");
		}
	
	/* package method */ void writeRDFXml(final XMLStreamWriter w) throws XMLStreamException {
		if(!isString()) {
			final String pfx = StringUtils.ifBlank(w.getPrefix(RDF.NS), RDF.pfx);	
			w.writeAttribute(pfx, RDF.NS, "dataType", this.getDatatypeURI());
			}
		else if(!StringUtils.isBlank(getLang())) {
			w.writeAttribute(
				XMLConstants.XML_NS_PREFIX,
				XMLConstants.XML_NS_URI,
				"lang",
				this.getLang()
				);
			}
		w.writeCharacters(this.getLiteralAsString());	
		}

	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	/* pacakage */ int _compareTo(final Literal o) {
		if(o==this) return 0;
		
		if(isNumber() && o.isNumber()) {
			if(isInteger() && o.isInteger()) {
				return asBigInteger().compareTo(o.asBigInteger());
				}
			return asBigDecimal().compareTo(o.asBigDecimal());
			}
		
		final Class<?> C1 = this.getValue().getClass();
		final Class<?> C2 = o.getValue().getClass();
		if(!C1.equals(C2)) {
			return C1.getName().compareTo(C2.getName());
			}
		final Comparable v1 = this.getValue();
		final Comparable v2 = o.getValue();
		int i =  v1.compareTo(v2);
		if(i==0 && C1.equals(String.class)) {
			i = getLang().compareTo(o.getLang());
			}
		return i;
		}
	
	@Override
	public int hashCode() {
		return Objects.hash(value,lang);
		}
	
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null || !(obj instanceof Literal))
			return false;
		final Literal other = (Literal) obj;
		if(isNumber() && other.isNumber()) {
			if(isInteger() && other.isInteger()) {
				return asBigInteger().equals(other.asBigInteger());
				}
			return asBigDecimal().equals(other.asBigDecimal());
			}
		
		return Objects.equals(getValue(), other.getValue()) && 
				(isString()? getLang().equals(other.getLang()) : true);
		}
	
	@Override
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		if(isString()) {
			sb.append(StringUtils.doubleQuote(this.getString()));
			if(hasLang()) {
				sb.append("(").append(this.lang).append(")");
				}
			}
		else
			{
			sb.append(this.value.toString()).append("@").append(getDatatypeURI());
			}
		return sb.toString();
		}
	}
