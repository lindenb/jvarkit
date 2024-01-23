package com.github.lindenb.jvarkit.tools.genome2svg;

import java.text.DecimalFormat;

import org.apache.commons.lang.StringUtils;
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
	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");

	
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
	
	public Element rect(Locatable locatable,double y,double height) {
		return rect(locatable.getStart(),locatable.getEnd(),y,height);
		}
	
	public Element rect(int chromStart1,int chromEnd1,double y,double height) {
		final Element rect = element("rect");
		double x1 = pos2pixel(trimpos(chromStart1));
		double x2 = pos2pixel(trimpos(chromEnd1+1));
		rect.setAttribute("x", format(x1));
		rect.setAttribute("y", format(y));
		rect.setAttribute("width", format(x2-x1));
		rect.setAttribute("height", format(height));
		return rect;
		}
	
	public Element anchor(Element child,final String url) {
		if(StringUtils.isBlank(url)) return child;
		final Element a = element("a");
		a.setAttribute("href", url);
		a.appendChild(child);
		return a;
		}
	
	public void appendStyle(String s) {
		if(StringUtils.isBlank(s)) return ;
		this.styleNode.appendChild(text(s+"\n"));
		}
	
	public String format(double v) {
		return this.decimalFormater.format(v);
		}
}
