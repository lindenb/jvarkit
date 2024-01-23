package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.text.DecimalFormat;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;


public abstract class Track {
	private static int ID_GENERATOR = 0;
	private final String id = "id"+(++ID_GENERATOR);
	protected String shortDesc = null;
	protected String longDesc = null;
	private String display = "";
	private String reference = null;
	private double featureHeight = 10;
	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	
	protected Track() {
		}
	
	public String getId() {
		return id;
		}
	
	public void setReference(String reference) {
		this.reference = reference;
		}
	public String getReference() {
		return reference;
		}

	public void setFeatureHeight(double featureHeight) {
		this.featureHeight = featureHeight;
		}
	public double getFeatureHeight() {
		return featureHeight;
		}
	
	
	public void setShortDesc(String shortDesc) {
		this.shortDesc = shortDesc;
	}
	
	public String getShortDesc() {
		return StringUtils.isBlank(this.shortDesc)?getClass().getName():this.shortDesc;
		}
	
	public void setLongDesc(String longDesc) {
		this.longDesc = longDesc;
		}
	public String getLongDesc() {
		return StringUtils.isBlank(this.longDesc)?getShortDesc():this.longDesc;
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
	
	/** convert double to string */
	protected String format(final double v)
		{
		return this.decimalFormater.format(v);
		}
	}
