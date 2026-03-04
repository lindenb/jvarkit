package com.github.lindenb.jvarkit.chart;

public class AbstractDataY {
	private double y;
	protected AbstractDataY(final double y) {
		this.y = y;
		}
	
	public double getY() {
		return y;
		}
	public void setY(double y) {
		this.y = y;
		}
	}
