package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.text.DecimalFormat;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class AbstractBean {
	private static int ID_GENERATOR = 0;
	private static final DecimalFormat decimalFormater = new DecimalFormat("##.##");

	private final String id = "id"+(++ID_GENERATOR);
	protected String shortDesc = null;
	protected String longDesc = null;

	protected AbstractBean() {
	}

	public String getId() {
		return id;
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
	
	/** convert double to string */
	protected String format(final double v)
		{
		return decimalFormater.format(v);
		}
	}
