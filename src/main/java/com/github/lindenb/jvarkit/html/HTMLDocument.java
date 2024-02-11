package com.github.lindenb.jvarkit.html;

import java.util.List;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.xml.DocumentWrapper;

public class HTMLDocument extends DocumentWrapper {
	public final Document document;
	public final Element htmlElement;
	public final Element headElement;
	public final Element titleElement;
	public final Element styleElement;
	public final Element scriptElement;
	public final Element bodyElement;
	public HTMLDocument() {
		this.document = DocumentWrapper.makedoc(true);
		this.htmlElement = element("html",Maps.of("lang", "en-US"));
		
		
		this.document.appendChild(htmlElement);
	
		this.headElement = element("head");
		this.htmlElement.appendChild(headElement);
		
		this.headElement.appendChild(element("meta",Maps.of("charset", "utf-8")));
		
		this.titleElement = element("title","untitled");
		this.headElement.appendChild(titleElement);
		
		this.styleElement = element("style");
		this.headElement.appendChild(styleElement);
		
		this.scriptElement = element("script");
		this.headElement.appendChild(scriptElement);

		this.bodyElement = element("body");
		this.htmlElement.appendChild(bodyElement);
	}
	
	public Node importSVG(SVGDocument dom) {
		return getDocument().importNode(dom.svgElement,true);
	}
	
	public Document getDocument() {
		return document;
		}
	
	private Element wrap(String tag,Node n) {
		final Element E = element(tag);
		E.appendChild(n);
		return E;
		}
	
	public Element bold(Node n) { return wrap("b",n);}
	public Element italic(Node n) { return wrap("i",n);}
	public Element th(Node n) { return wrap("th",n);}
	public Element td(Node n) { return wrap("td",n);}

	public Element anchor(Element wrapped,String url) {
		if(StringUtils.isBlank(url)) return wrapped;
		final Element a = element("a",Maps.of("href", url));
		a.appendChild(wrapped);
		return a;
		}

	public void setTitle(final String s) {
		removeAllChildren(this.titleElement).appendChild(text(s));
		}

	
	@Override
	public String getDefaultNamespace() {
		return "http://www.w3.org/1999/xhtml";
		}
	
	public class Table {
		Element table;
		Element thead;
		Element tbody;
		final int nCols;
		Table(List<Object> header) {
			this.table = element("table");
			this.thead = element("thead");
			this.table.appendChild(thead);
			this.tbody = element("tbody");
			this.table.appendChild(tbody);
			this.nCols = header.size();
			for(Object o:header) {
				final Element tr=element("tr");
				thead.appendChild(tr);
				for(Object t:header) {
					final Element th=element("th");
					tr.appendChild(th);
					th.appendChild(convert(t));
					}
				}
			}
		public Element get(int y) {
			int row=0;
			for(Node n = this.tbody.getFirstChild();n!=null ;n=n.getNextSibling(),row++)
				{
				if(row==y) return Element.class.cast(n);
				}
			}
		public Element get(int y,int x) {
			Element row=get(y);
			return null;
			}
		public void set(int y,int x, Object o) {
			Element r= get(y,x);
			removeAllChildren(r);
			r.appendChild(convert(o));
			}
		
		}
	
	
	
}
