package com.github.lindenb.jvarkit.tools.genome2svg.beans;


import java.awt.Color;
import java.util.Collections;
import java.util.List;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Log;


public abstract class Track extends AbstractBean {
	private String display = "";
	private double featureHeight = 10;
	private List<? extends PathBean> paths = Collections.emptyList();
	private Log logger = null;
	protected Track() {
		setLogger(htsjdk.samtools.util.Log.getInstance(getClass()));
		}

	private ReferenceBean reference = null;

	public void setReference(ReferenceBean reference) {
		this.reference = reference;
		}
	public ReferenceBean getReference() {
		return reference;
		}

	public Log getLogger() {
		return logger;
		}
	public void setLogger(Log logger) {
		this.logger = logger;
		}
	
	public void setPath(PathBean onePath) {
		if(onePath==null) {
			setPaths(Collections.emptyList());
			}
		else
			{
			setPaths(Collections.singletonList(onePath));
			}
		}

	public PathBean getPath() {
		switch(this.paths.size()) {
			case 0: return null;
			case 1: return this.paths.get(0);
			default: throw new IllegalStateException("Asking one path for but got "+this.paths.size());
			}
		}
	
	public void setPaths(List<? extends PathBean> onePath) {
		if(onePath==null) {
			this.paths =  Collections.emptyList();
			}
		else
			{
			this.paths = onePath;
			}
		}
	public List<? extends PathBean> getPaths() {
		return paths;
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
	

	protected String getColorPalette(float f) {
		if(f<0f || f>1f) return "none";
		return ColorUtils.toRGB( Color.getHSBColor(f, 0.85f, 1.0f));
		}
	
	protected String getUrlForLocatable(final Locatable loc) {
		return null;
		}
	
	}
