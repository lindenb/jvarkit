package com.github.lindenb.jvarkit.tools.genome2svg.beans;


import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;

import htsjdk.samtools.util.Locatable;


public abstract class Track extends AbstractBean {
	private String display = "";
	private double featureHeight = 10;
	
	
	protected Track() {
		}

	private ReferenceBean reference = null;

	public void setReference(ReferenceBean reference) {
		this.reference = reference;
		}
	public ReferenceBean getReference() {
		return reference;
		}


	public void setFeatureHeight(double featureHeight) {
		this.featureHeight = featureHeight;
		}
	public double getFeatureHeight() {
		return featureHeight;
		}

	
	public void setDisplay(String display) {
		this.display = display;
		}
	public String getDisplay() {
		return display;
	}
	public boolean isVisible() {
		return !getDisplay().equals("hidden");
		}
	abstract public void paint(SVGContext ctx);
	
	protected void insertTitle(Element root,SVGContext ctx) {
		String t = getShortDesc();
		if(StringUtils.isBlank(t)) return;
		int fontSize = 10;
		ctx.y += 2;
		Element txt = ctx.element("text");
		txt.setAttribute("x", format(ctx.image_width/2.0));
		txt.setAttribute("y", format(ctx.y+fontSize));
		txt.setAttribute("style", "stroke:none;fill:darkgray;stroke-width:1px;text-anchor:middle;font-size:"+fontSize+";");
		txt.appendChild(ctx.text(t));
		t = getLongDesc();
		if(!StringUtils.isBlank(t)) {
			txt.appendChild(ctx.title(t));
			}
		root.appendChild(txt);
		ctx.y += fontSize;
		ctx.y += 5;
		}
	

	protected String getUrlForLocatable(final Locatable loc) {
		return null;
		}
	
	
	
	}
