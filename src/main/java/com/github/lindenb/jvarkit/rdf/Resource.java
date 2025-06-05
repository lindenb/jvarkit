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

import java.net.URI;
import java.util.Objects;

import javax.xml.namespace.QName;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.RDFS;

public class Resource implements RDFNode , Comparable<Resource>{
	public static final Resource  RDF_type = new Resource(RDF.NS, "type");
	public static final Resource  RDFS_label = new Resource(RDFS.NS, "label");
	private static long ID_GENERATOR=System.currentTimeMillis();
	private final String namespaceURI;
	private final String localName;
	public Resource(String namespaceURI,String localName) {
		this.namespaceURI = namespaceURI;
		this.localName = localName;
		}
	
	public Resource(final String uri) {
		String[] tokens= splitURI(uri);
		this.namespaceURI= tokens[0];
		this.localName= tokens[1];
		}
	
	public Resource(QName qName) {
		this(qName.getNamespaceURI(),qName.getLocalPart());
		}
	public Resource() {
		this("blank:","_"+(++ID_GENERATOR));
		}
	
	
	public String getURI() {
		return this.namespaceURI+this.localName;
	}
	
	public String getLocalName() {
		return this.localName;
		}
	
	public String getNamespaceURI() {
		return this.namespaceURI;
		}
	
	private static String[] splitURI(final String uri) {
		if(StringUtils.isBlank(uri)) throw new IllegalArgumentException("blank uri");
		final String ns;
		final String lcl;
		int i = uri.lastIndexOf('#');
		if(i!=-1) {
			ns = uri.substring(0,i+1);
			lcl = uri.substring(i+1);
			} else
				{
				i = uri.lastIndexOf('/');
				if(i!=-1) {
					ns = uri.substring(0,i+1);
					lcl = uri.substring(i+1);
					} else
						{
						throw new IllegalArgumentException("cannot find '/' or '#' in uri");
						}
				}
		return new String[] {ns,lcl};
		}
	
	public static String namespaceOfURI(final String uri) {
		return splitURI(uri)[0];
		}
	
	public static String localNameOfURI(final String uri) {
		return splitURI(uri)[1];
		}
	
	public static boolean isURI(final String uri) {
		try {
			URI.create(uri);
			return true;
			}
		catch(Throwable err) {
			return false;
			}
		}
	
	@Override
	public boolean isResource() {
		return true;
		}
	@Override
	public boolean isLiteral() {
		return false;
		}
	
	@Override
	public int hashCode() {
		return Objects.hash(localName, namespaceURI);
		}
	
	
	@Override
	public int compareTo(Resource other) {
		int i= namespaceURI.compareTo(other.namespaceURI);
		if(i!=0) return i;
		return localName.compareTo(other.localName);
		}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Resource other = (Resource) obj;
		return namespaceURI.equals(other.namespaceURI) && 
				this.localName.equals(other.localName);
		}
	
	@Override
	public String toString() {
		return new StringBuilder("<").append(namespaceURI).append(localName).append(">").toString();
		}
	}
