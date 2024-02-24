package com.github.lindenb.jvarkit.html;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
		
		this.styleElement = element("style",
				getDefaultCsstyle().entrySet().stream().
					map(KV->KV.getKey()+" {"+KV.getValue()+"}").
					collect(Collectors.joining("\n")));
		this.headElement.appendChild(styleElement);
		
		this.scriptElement = element("script");
		this.headElement.appendChild(scriptElement);

		this.bodyElement = element("body");
		this.htmlElement.appendChild(bodyElement);
	}
	
	
	protected Map<String,String> getDefaultCsstyle() {
		final Map<String,String> h = new HashMap<>();
		h.put("table","border-collapse: collapse;");
		h.put("tbody tr:nth-child(odd)","background-color: #fff;");
		h.put("tbody tr:nth-child(even)","background-color: #eee;");
		h.put("caption","fbackground-color: #333;color: white;");
		h.put("thead","font-weight: bold;");
		h.put("td","border: 1px solid black");
		h.put("th","text-align:center;border: 1px solid black;");
		h.put("footer","text-align: center;  padding: 5px; background-color: #abbaba;color: #000;");
		return h;
		}
	
	public Node importSVG(SVGDocument dom) {
		return getDocument().importNode(dom.svgElement,true);
	}
	
	public Document getDocument() {
		return document;
		}
	
	
	public Element bold(Object content) { return element("b",content,null);}
	public Element italic(Object content) { return element("i",content,null);}
	public Element code(Object content) { return element("code",content,null);}
	public Element th(Object content) { return element("th",content,null);}
	public Element td(Object content) { return element("td",content,null);}

	public Node anchor(Object content,String url) {
		return element("a",content,Maps.of("href", url));
		}
	
	public Node anchor(String url) {
		return anchor(url,url);
		}

	public void setTitle(final String s) {
		removeAllChildren(this.titleElement).appendChild(text(s));
		}

	
	@Override
	public String getDefaultNamespace() {
		return "http://www.w3.org/1999/xhtml";
		}
	
	
	public class Table {
		public Element table;
		final Element thead;
		final Element tbody;
		public Element rowFoot=null;
		final int nCols;
		final Map<String,Integer> col2index= new HashMap<>();
		Table(List<Object> header) {
			this.table = element("table");
			this.thead = element("thead");
			this.table.appendChild(thead);
			this.tbody = element("tbody");
			this.table.appendChild(tbody);
			this.nCols = header.size();
			//
			final Element tr=element("tr");
			thead.appendChild(tr);
			for(Object t:header) {
				final Element th=element("th");
				tr.appendChild(th);
				Node lblNode = toNode(t);
				final String label = lblNode.getTextContent();
				col2index.put(label, col2index.size());
				th.appendChild(lblNode);
				}
				
			}
		
		
		public int getColumnCount() {
			return nCols;
			}
		
		public int getRowCount() {
			int row=0;
			for(Node n = this.tbody.getFirstChild();n!=null ;n=n.getNextSibling())
				{
				row++;
				}
			return row;
			}
		
		public Element appendRow() {
			return get(getRowCount());
			}
		
		private Element createRow() {
			final Element row = element("tr");
			for(int i=0;i< getColumnCount();++i) {
				final Element td=element("td");
				row.appendChild(td);
				}
			return row;
			}
		
		public Element get(int y) {
			int row=0;
			for(Node n = this.tbody.getFirstChild();n!=null ;n=n.getNextSibling())
				{
				if(row==y) return Element.class.cast(n);
				row++;
				}
			Element last=null;
			while(row<= y) {
				last =createRow();
				this.tbody.appendChild(last);
				row++;
				}
			return last;
			}
		public Element get(int y,int x) {
			int c=0;
			Element row=get(y);
			for(Node n = row.getFirstChild();n!=null ;n=n.getNextSibling())
				{
				if(c==x) return Element.class.cast(n);
				c++;
				}
			return null;
			}
		public Element get(int y,String label) {
			return get(y,this.col2index.getOrDefault(label,-1));
			}
		public void set(int y,int x, Object o) {
			Element r= get(y,x);
			removeAllChildren(r);
			r.appendChild(toNode(o));
			}
		public void set(int y,String label, Object o) {
			set(y,this.col2index.getOrDefault(label,-1),o);
			}
		public void setFooter(int x, Object o) {
			if(this.rowFoot==null) {
				this.rowFoot = createRow();
				Element e= element("tfoot");
				e.appendChild(rowFoot);
				this.table.appendChild(e);
				}
			int c=0;
			for(Node n = rowFoot.getFirstChild();n!=null ;n=n.getNextSibling())
				{
				if(c==x) {
					removeAllChildren(n);
					n.appendChild(toNode(o));
					break;
					}
				c++;
				}
			}
		}
	
	public Element div() {
		return element("div");
		}
	
	public Element hr() {
		return element("hr");
		}
	public Element paragraph() {
		return element("p");
		}

	
	public Element header(int level) {
		return element("h"+level);
		}
	
	public Element h1() { return header(1); }
	public Element h2() { return header(2); }
	public Element h3() { return header(3); }
	public Element h4() { return header(4); }
	
	public Element span(Object content,Map<String,Object> atts) {
		Element  b = element("span",atts);
		if(content!=null) b.appendChild(toNode(content));
		return b;
		}

	
	public Element button(Object content,Map<String,Object> atts) {
		Element  b = element("button",atts);
		if(content!=null) b.appendChild(toNode(content));
		return b;
		}
	
	public Element textarea(String content,Map<String,Object> atts) {
		return  element("textarea",content,atts);
		}
	
	public Table createTable(List<Object> header) {
		return new Table(header);
	}
	
	
}
