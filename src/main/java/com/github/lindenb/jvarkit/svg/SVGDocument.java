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
package com.github.lindenb.jvarkit.svg;

import java.awt.geom.Point2D;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.lang.primitive.DoubleArray;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.xml.DocumentWrapper;

public class SVGDocument extends DocumentWrapper {
public final Document document;
public final Element svgElement;
public final Element styleElement;
public final Element titleElement;
public final Element metaElement;
public final Element rootElement;
public final Element defsElement;
public final Element scriptElement;
public final Element rdfRoot;
private final Map<String, String> _style2class = new HashMap<>();
/**
 * SVGDocument
 * @param document
 */
public SVGDocument(final Document document) {
	this.document =  Objects.requireNonNull(document);
	this.svgElement = element( "svg",Maps.of("version","1.1"));
	this.document.appendChild(this.svgElement);
	this.svgElement.appendChild(document.createComment("jvarkit.version: "+JVarkitVersion.getInstance().getGitHash()));

	
	this.titleElement = element("title","untitled");
	this.svgElement.appendChild(this.titleElement);

	this.metaElement = element("metadata");
	this.svgElement.appendChild(this.metaElement);
	
	
	this.rdfRoot = document.createElementNS(RDF.NS, "rdf:RDF");
	this.metaElement.appendChild(this.rdfRoot);

	
	this.styleElement = element("style",
			"svg {stroke:black;fill:none;stroke-width:1px;}\n"+
					"rect {stroke:black;fill:none;}\n"+
					"polygon {stroke:black;fill:none;}\n"+
					"path {stroke:black;fill:none;}\n"+
					"polyline {stroke:black;fill:none;}\n"+
					"text {stroke:none;fill:black;}\n"
			);
	this.svgElement.appendChild(this.styleElement);
	
	
	this.defsElement = element("defs");
	this.svgElement.appendChild(this.defsElement);

	this.scriptElement = element("script");
	this.svgElement.appendChild(this.scriptElement);

	
	this.rootElement = element("g");
	this.svgElement.appendChild(this.rootElement);
	
	}

public SVGDocument() {
	this(DocumentWrapper.makedoc(true));
	}


@SuppressWarnings("serial")
private static final Set<String> LOOKS_LIKE_STYLE = new HashSet<String>() {{{
	add("stroke");
	add("fill");
	add("text-anchor");
	add("text-align");
	add("opacity");
	add("font-size");
	add("stroke-width");
	add("opacity");
	}}};

@Override
protected Map<String, Object> updateAttributes(Map<String, Object> atts) {
	if(atts==null || atts.containsKey("style") || atts.containsKey("class")) return atts;
	Map<String,String> style= null;
	for(final String key:atts.keySet()) {
		if(LOOKS_LIKE_STYLE.contains(key)) {
			if(style==null) style = new HashMap<>();
			style.put(key, convertObjectToString(atts.get(key)));
			}
		}
	if(style==null) {
		return atts;
		}
	//copy if unmodifiable
	atts = new HashMap<>(atts);
	for(final String key:style.keySet()) {
		atts.remove(key);
		}
	atts.put("style", style.entrySet().stream().map(KV->KV.getKey()+":"+KV.getValue()).collect(Collectors.joining(";")));
	return atts;
	}

public void appendStyle(final String s) {
	this.styleElement.appendChild(text(s+"\n"));
	}

public void setTitle(final String s) {
	removeAllChildren(this.titleElement).appendChild(text(s));
	}
public void setWidth(final double w) {
	this.svgElement.setAttribute("width", format(w));
	}
public double getWidth() {
	String s= this.svgElement.getAttribute("width");
	return StringUtils.isBlank(s)?0:(int)Double.parseDouble(s);
	}


public void setHeight(final double h) {
	this.svgElement.setAttribute("height", format(h));
	}

public double getHeight() {
	String s= this.svgElement.getAttribute("height");
	return StringUtils.isBlank(s)?0:(int)Double.parseDouble(s);
	}

public Document getDocument() {
	return document;
	}
@Override
public String getDefaultNamespace() {
	return SVG.NS;
	}

/** if given style is not delcared in the style element, add it to style and return the class-id */
public String style2class(final String style) {
	return style2class(null,style);
	}


/** if given style is not delcared in the style element, add it to style and return the className */
public String style2class(final String className,final String style) {
	String id = this._style2class.get(style);
	if(id!=null) return id;
	id = StringUtils.ifBlank(className, "c"+nextId());
	appendStyle("."+id+"{"+style+"}");
	this._style2class.put(style, id);
	return id;
	}


public Element line(double x1,double y1,double x2,double y2,Map<String,Object> atts) {
	final Map<String,Object> m =Maps.of("x1",x1,"y1",y1,"x2",x2,"y2",y2);
	if(atts!=null) m.putAll(atts);
	return  this.element("line",null,m);
	}

public Element line(double x1,double y1,double x2,double y2) {
	return this.line(x1, y1, x2, y2,Collections.emptyMap());
	}


public Element rect(double x,double y,double width,double height,Map<String,Object> atts) {
	final Map<String,Object> m = Maps.of("x",x,"y",y,"width",width,"height",height);
	if(atts!=null) m.putAll(atts);
	return  this.element("rect",null,m);
	}

public Element rect(double x,double y,double width,double height) {
	return this.rect(x, y, width, height,Collections.emptyMap());
	}


public Element text(double x,double y,Object str,Map<String,Object> atts) {
	final Map<String,Object> m =  Maps.of("x", x, "y", y);
	if(atts!=null) m.putAll(atts);
	return this.element("text",str,m);
	}

public Element text(double x,double y,Object str) {
	return this.text(x, y, str,Collections.emptyMap());
	}

public Element circle(double cx,double cy,double r, final Map<String,Object> atts) {
	final Map<String,Object> m = Maps.of("cx", cx, "cy", cy, "r", r);
	if(atts!=null) m.putAll(atts);
	return  this.element("circle",null,m);
	}

public Element circle(double cx,double cy,double r) {
	return this.circle(cx, cy, r,Collections.emptyMap());
	}

public Element use(double x,double y,String id, final Map<String,Object> atts) {
	final Map<String,Object> m = Maps.of("x", x, "y", y, "href", (id.startsWith("#")?id:"#"+id));
	if(atts!=null) m.putAll(atts);
	return  this.element("use",null,m);
	}

public Element use(double x,double y,String id) {
	return this.use(x, y, id,Collections.emptyMap());
	}


public Element group(Map<String,Object> atts) {
	return this.element("g",null,atts);
	}

public Element clipPath() {
	return this.element("clipPath");
	}

public Element polyline(final Map<String,Object> atts) {
	return this.element("polyline",null,atts);
	}
public Element polyline() {
	return polyline(Collections.emptyMap());
	}

public String toString(final List<? extends Point2D> points) {
	final StringBuilder sb = new StringBuilder();
	for(int i=0;i< points.size();i++) {
		if(i>0) sb.append(" ");
		sb.append(format(points.get(i).getX()));
		sb.append(",");
		sb.append(format(points.get(i).getY()));
		}
	return sb.toString();
	}

public Element group(Node...items) {
	final Element g = group();
	for(Node n:items) {
		g.appendChild(n);
		}
	return g;
	}


public Element group() {
	return group(Collections.emptyMap());
	}

/** return a css selector for translate */
public String translate(double dx,double dy) {
	return "translate("+format(dx)+","+format(dy)+")";
}

public class Shape implements Supplier<Element> {
	private final StringBuilder sb = new StringBuilder();
	private Shape append(char symbol,double...array) {
		if(sb.length()>0) sb.append(" ");
		sb.append(symbol);
		for(int i=0;i< array.length;++i) {
			sb.append(i==0?" ":",");
			sb.append(format(array[i]));
			}
		return this;
		}
	public Shape moveTo(double x,double y) {
		return append('M',x,y);
		}
	public Shape moveToR(double x,double y) {
		return append('m',x,y);
		}
	public Shape lineTo(double x,double y) {
		return append('L',x,y);
		}
	public Shape lineToR(double x,double y) {
		return append('l',x,y);
		}
	
	public Shape horizonal(double h) {
		return append('H',h);
		}
	public Shape horizonalR(double h) {
		return append('h',h);
		}
	public Shape vertical(double v) {
		return append('V',v);
		}
	public Shape verticalR(double v) {
		return append('v',v);
		}
	
	public Shape closePath() {
		return append('z');
		}
	
	public Shape bezier(double x1,double y1, double x2,double y2,double x,double y) {
		return append('C',x1,y1,x2,y2,x,y);
		}

	public Shape bezierR(double x1,double y1, double x2,double y2,double dx,double dy) {
		return append('c',x1,y1,x2,y2,dx,dy);
		}

	
	public Shape quadratic(double x1,double y1, double x,double y) {
		return append('Q',x1,y1,x,y);
		}

	public Shape quadraticR(double dx1,double dy1, double dx,double dy) {
		return append('q',dx1,dy1,dx,dy);
		}

	public Shape arc(double rx,double ry, double r_axis_rotation,boolean large_arc_flag,boolean sweep_flag,double x,double y) {
		return append('A',rx,ry,r_axis_rotation,(large_arc_flag?0:1),(sweep_flag?0:1),x,y);
		}

	public Shape arcR(double rx,double ry, double r_axis_rotation,boolean large_arc_flag,boolean sweep_flag,double x,double y) {
		return append('a',rx,ry,r_axis_rotation,(large_arc_flag?0:1),(sweep_flag?0:1),x,y);
		}

	
	public Element make(Map<String,Object> atts) {
		final Element e= SVGDocument.this.element("path",atts);
		e.setAttribute("d", sb.toString());
		return e;
		}
	
	public Element make() {
		return this.make(Collections.emptyMap());
		}

	@Override
	public Element get() {
		return this.make();
		}

	}

public class Polygon implements Supplier<Element> {
	private final DoubleArray xarray = new DoubleArray();
	private final DoubleArray yarray = new DoubleArray();
	private final String tagName;
	Polygon(final String tagName) {
		this.tagName = tagName;
		}
	public Polygon lineTo(double x,double y) {
		this.xarray.add(x);
		this.yarray.add(y);
		return this;
		}
	public Polygon horizontal(double h) {
		return lineTo(h,getLastY());
		}
	public Polygon vertical(double v) {
		return lineTo(getLastX(),v);
		}

	public Polygon horizontalRel(double h) {
		return horizontal(getLastX()+h);
		}
	public Polygon verticalRel(double v) {
		return vertical(getLastY()+v);
		}
	public double getLastX() {
		if(this.size()==0) throw new IllegalStateException("polygon is empty");
		return X(size()-1);
		}
	public double getLastY() {
		if(this.size()==0) throw new IllegalStateException("polygon is empty");
		return Y(size()-1);
		}

	private int size() { return xarray.size();}
	private double X(int i) { return xarray.get(i);}
	private double Y(int i) { return yarray.get(i);}
	
	
	private void simplify() {
		int i=0;
		while(i < size()) {
			if(i+1 < size() &&
			   X(i) == X(i+1) &&
			   Y(i) == Y(i+1))
				{
				xarray.remove(i+1);
				yarray.remove(i+1);
				continue;
				}
			else if(i+2 < size() && Y(i) == Y(i+1) && Y(i+1)==Y(i+2))
				{
				xarray.remove(i+1);
				yarray.remove(i+1);
				continue;
				}
			else if(i+2 < size() && X(i) == X(i+1) && X(i+1)==X(i+2))
				{
				xarray.remove(i+1);
				yarray.remove(i+1);
				continue;
				}

			else
				{
				i++;
				}
			}
		}
	
	
	
	public Element make(Map<String,Object> atts) {
		final StringBuilder sb = new StringBuilder();
		simplify();
		for(int i=0;i< this.size();i++) {
			if(i>0) sb.append(" ");
			sb.append(format(X(i)));
			sb.append(",");
			sb.append(format(Y(i)));
			}
		final Element e= SVGDocument.this.element(this.tagName,atts);
		e.setAttribute("points", sb.toString());
		return e;
		}
	public Element make() {
		return this.make(Collections.emptyMap());
		}

	@Override
	public Element get() {
		return this.make();
		}
	}


public Polygon createPolygon() {
	return new Polygon("polygon");
	}
public Polygon createPolyline() {
	return new Polygon("polyline");
	}

public Shape createShape() {
	return new Shape();
	}

public Element anchor(final Element wrapped,final String url) {
	if(StringUtils.isBlank(url)) return wrapped;
	final Node parent = wrapped.getParentNode();
	final Element a = element("a",Maps.of("href", url,"target","_blank"));
	if(parent!=null) {
		parent.replaceChild(a, wrapped);
		}
	a.appendChild(wrapped);
	return a;
	}

public Element setTitle(final Element root,final String  s) {
	if(root==null || StringUtils.isBlank(s)) return root;
	final Element e= element("title",s);
	root.appendChild(e);
	return root;
	}

public String createVerticalLinearGradient(String color1,String color2) {
	Node n= this.makeNodeBuilder("linearGradient").
		attribute("x1","50%").
		attribute("x2","50%").
		attribute("y1","0%").
		attribute("y2","100%").
		startElement("stop").
			attribute("offset","0%").
			attribute("style","stop-color:"+color1+";stop-opacity:1;").
		endElement().
		startElement("stop").
			attribute("offset","50%").
			attribute("style","stop-color:"+color2+";stop-opacity:1;").
		endElement().
		startElement("stop").
			attribute("offset","100%").
			attribute("style","stop-color:"+color1+";stop-opacity:1;").
		endElement().
		make();
	return insertDefElement(Element.class.cast(n));
	}
/** insert element in <defs> , check if very same element exists add new @id if missing, return @id */
public String insertDefElement(final Element e) {
	for(Node c1=this.defsElement.getFirstChild();c1!=null;c1=c1.getNextSibling()) {
		if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
		Element e1 = Element.class.cast(c1);
		if(!e.hasAttribute("id")) continue;
		final String id = e.getAttribute("id");
		if(StringUtils.isBlank(id)) continue;
		if(same(e1,e,0)) {
			return id;
			}
		}
	
	if(!e.hasAttribute("id")) e.setAttribute("id", nextId());
	this.defsElement.appendChild(e);
	return e.getAttribute("id");
	}

private boolean same(Element root, Element other,int depth) {
	if(!root.getNodeName().equals(other.getNamespaceURI())) return false;
	NamedNodeMap nm1 = root.getAttributes();
	NamedNodeMap nm2 = other.getAttributes();
	for(int i=0;i< nm1.getLength();i++) {
		Attr att1 = (Attr)nm1.item(i);
		if(depth==0 && att1.getNodeName().equals("id")) continue;
		Attr att2=  (Attr)nm2.getNamedItem(att1.getName());
		if(att2==null) return false;
		if(!att1.getNodeValue().equals(att2.getNodeValue())) return false;
		}
	for(int i=0;i< nm2.getLength();i++) {
		Attr att1 = (Attr)nm2.item(i);
		if(depth==0 && att1.getNodeName().equals("id")) continue;
		Attr att2=  (Attr)nm1.getNamedItem(att1.getName());
		if(att2==null) return false;
		if(!att1.getNodeValue().equals(att2.getNodeValue())) return false;
		}
	final NodeList L1 = root.getChildNodes();
	final NodeList L2 = other.getChildNodes();
	if(L1.getLength()!=L2.getLength()) return false;
	for(int i=0;i<L1.getLength();i++) {
		Node c1= L1.item(i);
		Node c2= L2.item(i);
		if(c1.getNodeType()!=c2.getNodeType()) return false;
		if(c1.getNodeType()==Node.COMMENT_NODE) continue;
		else if(c1.getNodeType()==Node.ELEMENT_NODE) {
			if(!same(Element.class.cast(c1),Element.class.cast(c2),depth+1)) return false;
			}
		else if(c1.getNodeType()==Node.TEXT_NODE) {
			if(!Text.class.cast(c1).getTextContent().equals(Text.class.cast(c2).getTextContent())) return false;
			}
		else
			{
			return false;
			}
		}
	return true;
	}
}
