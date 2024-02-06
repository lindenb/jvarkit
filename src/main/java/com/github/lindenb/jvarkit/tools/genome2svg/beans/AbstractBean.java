package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.text.DecimalFormat;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class AbstractBean {
	private static int ID_GENERATOR = 0;
	private static final DecimalFormat decimalFormater = new DecimalFormat("##.##");

	private final String id = "id"+(++ID_GENERATOR);
	protected String shortDesc = null;
	protected String longDesc = null;
	protected String stroke=null;
	protected String fill=null;
	protected AbstractBean() {
	}

	public String getId() {
		return id;
		}
	
	
	
	public void setShortDesc(String shortDesc) {
		this.shortDesc = shortDesc;
	}
	
	public String getShortDesc() {
		return  this.shortDesc;
		}
	
	public void setLongDesc(String longDesc) {
		this.longDesc = longDesc;
		}
	public String getLongDesc() {
		return this.longDesc;
		}
	
	/** convert double to string */
	protected String format(final double v)
		{
		return decimalFormater.format(v);
		}
	
	public void setFill(String fill) {
		this.fill = fill;
		}
	public String getFill() {
		return fill;
		}
	public void setStroke(String stroke) {
		this.stroke = stroke;
		}
	
	public String getStroke() {
		return stroke;
		}
	
	@Override
	public String toString() {
		return String.valueOf(getClass())+" id:"+getId();
		}
	}
