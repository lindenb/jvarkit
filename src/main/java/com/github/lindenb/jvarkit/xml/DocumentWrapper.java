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
package com.github.lindenb.jvarkit.xml;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.StringWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Result;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.TransformerFactoryConfigurationError;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Attr;
import org.w3c.dom.CDATASection;
import org.w3c.dom.Comment;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.Text;

import com.github.lindenb.jvarkit.lang.StringUtils;


/**
 * DocumentWrapper
 * @author lindenb
 *
 */
public abstract class DocumentWrapper  {
	public abstract Document getDocument();
	private static int ID_GENERATOR = 0;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private Function<Object,String> stringConverter = new Function<Object,String>() {
		@Override
		public String apply(Object s) {
			if(s==null) return "";
			if(s instanceof Double) return format(Double.class.cast(s).doubleValue());
			if(s instanceof Float) return format(Float.class.cast(s).doubleValue());
			return String.valueOf(s);
			}
		};
	/** utiliy class to build a Dom tree */
	public class NodeBuilder {
		private NodeBuilder parent;
		private Node root;
		
		NodeBuilder() {
			this.root = DocumentWrapper.this.fragment();
			}
		NodeBuilder(String name,Map<String,Object> hash) {
			this.root = DocumentWrapper.this.element(name,hash);
			}
		
		public NodeBuilder startElement(String name,Map<String,Object> hash) {
			NodeBuilder b=new  NodeBuilder(name,hash);
			b.parent = this;
			this.root.appendChild(b.root);
			return b;
			}
		
		public NodeBuilder startElement(String name) {
			return startElement(name,Collections.emptyMap());
			}
		
		public NodeBuilder text(Object c) {
			if(c!=null) this.root.appendChild(DocumentWrapper.this.text(c));
			return this;
			}
		public NodeBuilder comment(Object c) {
			if(c!=null) this.root.appendChild(DocumentWrapper.this.comment(c)); 
			return this;
			}
		public NodeBuilder attribute(String key,Object v) {
			if(v!=null) {
				Element.class.cast(this.root).setAttribute(key,DocumentWrapper.this.convertObjectToString(v));
				}
			return this;
			}
		public NodeBuilder endElement() {
			return this.parent;
			}
		public Node make() {
			if(this.parent!=null) return parent.make();
			return this.root;
			}
		public Element makeElement() {
			final Node x = make();
			if(x.getNodeType()!=Node.ELEMENT_NODE) throw new IllegalStateException();
			return Element.class.cast(x);
			}
		}
	
	
	public NodeBuilder makeNodeBuilder() { return new NodeBuilder();}
	public NodeBuilder makeNodeBuilder(String name,Map<String,Object> att) { return new NodeBuilder(name,att);}
	public NodeBuilder makeNodeBuilder(String name) { return new NodeBuilder(name,Collections.emptyMap());}
	
	
	public String getDefaultNamespace() {
		return XMLConstants.NULL_NS_URI;
		}
	
	
	public String format(double v) {
		if((long)v==v) return String.valueOf((long)v);
		return this.decimalFormater.format(v);
		}
	
	public Function<Object,String> getStringConverter() {
		return this.stringConverter;
		}
	public void setStringConverter(Function<Object, String> stringConverter) {
		this.stringConverter = Objects.requireNonNull(stringConverter);
		}
	
	
	public String convertObjectToString(final Object s) {
		return getStringConverter().apply(s);
		}
	
	
	
	public Text text(final Object s) {
		return getDocument().createTextNode(convertObjectToString(s));
		}
	
	public Comment comment(final Object s) {
		return getDocument().createComment(convertObjectToString(s));
		}
	public CDATASection cdataSection(final Object s) {
		return getDocument().createCDATASection(convertObjectToString(s));
		}
	
	/** create a document fragment */
	public DocumentFragment fragment() {
		return getDocument().createDocumentFragment();
		}
	
	public Attr attribute(final String key) {
		return getDocument().createAttribute(key);
		}

	
	public Attr attribute(final String key,final Object v) {
		final Attr att= attribute(key);
		if(v!=null) att.setNodeValue(convertObjectToString(v));
		return att;
		}

	/** give a chance to change the map before it's used as a source of attribute. For example,
	 * in SVG change {stroke=x,fill=y} to {style=stroke:x;fill:y} 
	 * @param atts may be null
	 * @return
	 */
	protected Map<String,Object> updateAttributes( Map<String, Object> atts) {
		return atts;
		}
	
	public Element element(final String tag,Object content, Map<String, Object> atts) {
		final String ns = getDefaultNamespace();
		final Element e;
		if(StringUtils.isBlank(ns) || XMLConstants.NULL_NS_URI.equals(ns)) {
			e = getDocument().createElement(tag);
		} else {
			e = getDocument().createElementNS(ns, tag);
			}
		
		/* give a chance to change the attribute for example to group SVG stroke, fill, etc under 'style' */
		atts = updateAttributes(atts);
		
		if(atts!=null && !atts.isEmpty()) {
			for(final String key: atts.keySet()) {
				final Object v = atts.get(key);
				if(v==null) continue;
				e.setAttribute(key, convertObjectToString(v));
				}
			}
		
		if(content!=null) {
			e.appendChild(this.toNode(content));
			}
		return e;
		}
	public Element element(final String tag,CharSequence /* not object to avoid ambiguity */ textContent) {
		return this.element(tag, textContent, Collections.emptyMap());
	}
	
	public Element element(final String tag, final Map<String, Object> atts) {
		return this.element(tag, null, atts);
	}
	public Element element(final String tag) {
		return this.element(tag, null,  Collections.emptyMap());
		}
	
	protected <E extends Node> E setTextContent(final E n,Object o) {
		removeAllChildren(n).appendChild(text(o));
		return n;
		}
	
	protected <E extends Node> E removeAllChildren(final E n) {
		if(n==null) return n;
		while(n.hasChildNodes()) {
			n.removeChild(n.getFirstChild());
			}
		return n;
		}
	
	protected void saveTo(final Result out) throws IOException {
		try {
			TransformerFactory.newInstance().newTransformer().
				transform(new DOMSource(getDocument()),out);
		} catch (TransformerException | TransformerFactoryConfigurationError e) {
			throw new IOException(e);
			}
		}

	public void saveTo(final Writer out) throws IOException {
		saveTo(new StreamResult(out));
		}
	
	public void saveTo(final OutputStream out) throws IOException {
		saveTo(new StreamResult(out));
		}

	
	public void saveTo(final File out) throws IOException {
		saveTo(new StreamResult(out));
		}
	public void saveTo(final Path out) throws IOException {
		saveTo(out.toFile());
		}
	
	public void saveToFileOrStdout(final Path pathOrNull) throws IOException {
		if(pathOrNull!=null) {
			saveTo(pathOrNull);
			}
		else
			{
			saveTo(System.out);
			}
		}

	
	@Override
	public String toString() {
		try {
			final StringWriter sw = new StringWriter();
			saveTo(sw);
			return sw.toString();
		} catch (IOException e) {
			return "<?xml?>";
		}
		}
	
	public String nextId() {
		return "n"+(++ID_GENERATOR);
	}
	/** return o is it a org.w3c.dom.Node, DocumentFragment is List of array, otherwise, convert it to text */
	public Node toNode(Object o)  {
		if(o==null) return text("");
		if(o instanceof Node) {
			return Node.class.cast(o);
			}
		else if(o instanceof Collection) {
			final Collection<?> L =Collection.class.cast(o);
			final DocumentFragment df = this.fragment();
			for(Object c: L) {
				df.appendChild(toNode(c));
				}
			return df;
			}
		else if(o.getClass().isArray()) {
			return toNode(Arrays.asList((Object[])o));
			}
		else
			{
			return text(o);
			}
		}
	
	/** collect a stream of nodes and put them in a list, inserting sep as a text separator */
	public Collector<Node,?,DocumentFragment> joining(String sep) {
		return Collectors.collectingAndThen(
				Collectors.toList(),
				(L)->{
					final DocumentFragment f=fragment();
					for(int i=0;i< L.size();i++) {
						if(i>0 && sep!=null && !sep.isEmpty()) f.appendChild(text(sep));
						f.appendChild(L.get(i));
						}
					return f;
					}
				);
			}
	
	/** collect a stream of nodes and put them in a list, inserting sep as a text separator */
	public Collector<Node,?,DocumentFragment> joining() {
		return joining(null);
		}
	
	protected static Document makedoc(boolean namespaceAware) {
		try {
			final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			dbf.setNamespaceAware(namespaceAware);
			final DocumentBuilder db = dbf.newDocumentBuilder();
			return  db.newDocument();
		} catch (ParserConfigurationException e) {
			throw new RuntimeException(e);
			}
		}


}
