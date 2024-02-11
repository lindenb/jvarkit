/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
import java.util.Collections;
import java.util.Map;

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
	
	public String getDefaultNamespace() {
		return null;
		}
	
	public String convertObjectToString(final Object s) {
		return s==null?"":String.valueOf(s);
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
	
	public Attr attribute(final String key) {
		return getDocument().createAttribute(key);
		}

	
	public Attr attribute(final String key,final Object v) {
		final Attr att= attribute(key);
		if(v!=null) att.setNodeValue(convertObjectToString(v));
		return att;
		}

	
	public Element element(final String tag,Object textContent,final Map<String, Object> atts) {
		final String ns = getDefaultNamespace();
		final Element e;
		if(StringUtils.isBlank(ns)) {
			e = getDocument().createElement(tag);
		} else {
			e = getDocument().createElementNS(ns, tag);
			}
		
		if(atts!=null && !atts.isEmpty()) {
			for(final String key: atts.keySet()) {
				final Object v = atts.get(key);
				if(v==null) continue;
				e.setAttribute(key, convertObjectToString(v));
				}
			}
		
		if(textContent!=null) {
			e.appendChild(this.text(textContent));
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
	
	public Node toNode(Object o)  {
		if(o==null) return text("");
		if(o instanceof Node) {
			return Node.class.cast(o);
			}
		else
			{
			return text(o);
			}
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
