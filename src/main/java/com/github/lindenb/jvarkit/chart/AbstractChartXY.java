package com.github.lindenb.jvarkit.chart;

public abstract class AbstractChartXY extends Chart {
	private final Double[] y_limits = new Double[] {null,null};
	private String yAxisLabel="";
	private final Double[] x_limits = new Double[] {null,null};
	private String xAxisLabel="";

	public AbstractChartXY setYLimits(final Double m,Double M) {
		this.y_limits[0] = m;
		this.y_limits[1] = M;
		return this;
		}
	
	 public void setYAxisLabel(String yAxisLabel) {
			this.yAxisLabel = yAxisLabel;
		 	}
	
	 public String getYAxisLabel() {
		return yAxisLabel;
	 	}
	 

	
	public AbstractChartXY setXLimits(final Double m,Double M) {
		this.x_limits[0] = m;
		this.x_limits[1] = M;
		return this;
		}
	 public void setXAxisLabel(String xAxisLabel) {
		this.xAxisLabel = xAxisLabel;
	 	}
	 public String getXAxisLabel() {
		return xAxisLabel;
	 	}
}
