package com.github.lindenb.jvarkit.chart;

public abstract class Chart {
	public static abstract class AbstractData<X,Y>
		{
		private X x;
		private Y y;
		public AbstractData() {
			this(null,null);
			}
		public AbstractData(final X x,final Y y) {
			this.x = x;
			this.y = y;
			}
		public X getX() {
			return x;
			}
		public Y getY() {
			return y;
			}
		@Override
		public String toString() {
			return String.valueOf(this.x)+"/"+String.valueOf(this.y);
			}
		}
	
	private String legend = "";
	private boolean legendVisible = true;
	private String title = "";
	private String style = "";
	private boolean animated = false;
	
	
	public void setAnimated(boolean animated) {
		this.animated = animated;
	}
	
	public boolean isAnimated() {
		return animated;
	}
	
	public void setLegend(String legend) {
		this.legend = legend;
		}
	public String getLegend() {
		return legend;
		}
	
	public void setLegendVisible(boolean legendVisible) {
		this.legendVisible = legendVisible;
	}
	public boolean isLegendVisible() {
		return legendVisible;
	}
	
	public String getTitle() {
		return title;
		}
	public void setTitle(String title) {
		this.title = title;
		}
	public String getStyle() {
		return style;
	}
	public void setStyle(String style) {
		this.style = style;
	}
	public abstract void update();
}
