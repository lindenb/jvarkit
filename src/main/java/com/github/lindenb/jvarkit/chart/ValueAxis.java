package com.github.lindenb.jvarkit.chart;

import com.github.lindenb.jvarkit.chart.Axis.Type;

public abstract class ValueAxis<T extends Number> extends Axis<T> {
	private Double lowerBound = null;
	private Double upperBound = null;
	private Double ticks = null;
	
	protected ValueAxis(final String label) {
		this(label,null,null,null);
	}
	protected ValueAxis(final String label,Double lowerBound,Double upperBound,Double ticks) {
		super(label);
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.ticks = ticks;
	}
	protected ValueAxis(Double lowerBound,Double upperBound,Double ticks) {
		this("",lowerBound,upperBound,ticks);
	}
	
	public void setLowerBound(Double lowerBound) {
		this.lowerBound = lowerBound;
		}
	public Double getLowerBound() {
		return lowerBound;
	}
	public void setUpperBound(Double upperBound) {
		this.upperBound = upperBound;
	}
	public Double getUpperBound() {
		return upperBound;
	}
	
	@Override
	public Type getType() {
		return Type.NUMBER;
		}
}
