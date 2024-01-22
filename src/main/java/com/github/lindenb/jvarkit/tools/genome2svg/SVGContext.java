package com.github.lindenb.jvarkit.tools.genome2svg;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.util.Locatable;

public class SVGContext {
	public double y=0;
	public int left_margin = 100;
	public Document svgDom;
	public Element defsNode;
	public Element styleNode;
	public Element tracksNode;
	public Locatable loc;
	public double image_width;
	
	public Element element(final String name) {
		return svgDom.createElementNS(SVG.NS,name);
	}
	public Element element(final String name,final String content) {
		Element e = element(name);
		e.appendChild(text(content));
		return e;
		}
	
	public double pos2pixel(double pos1) {
		return ((pos1 - loc.getStart())/(double)loc.getLengthOnReference()) * image_width;
	}
	
	public double pixel2genomic(double pixX) {
		 return loc.getStart()+(pixX/image_width)*loc.getLengthOnReference();
	}
	
	public double trimpos(double pos1) {
		return Math.max(loc.getStart(), Math.min(loc.getEnd(), pos1));
	}
	public Text text(final String s) {
		return svgDom.createTextNode(s);
	}
	public Element title(final String s) {
		return element("title",s);
	}
}
