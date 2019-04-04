package com.github.lindenb.jvarkit.chart;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class PieChart extends Chart {
	
	public static class Data
		extends AbstractData<String,Number>
		{
		public Data(final String x,final double y) {
			super(x,y);
			if(y<0) throw new IllegalArgumentException("y Cannot be negative:"+y);
			}
		public String getName() {
			return getX();
			}
		public double getPieValue() {
			return getY().doubleValue();
			}
		}
	private final List<Data> data = new ArrayList<>();
	
	public PieChart(final List<Data> data) {
		this.data.addAll(data);
		}
	public PieChart() {
		this(Collections.emptyList());
		}
	
	public boolean isEmpty() {
		return getData().isEmpty();
	}
	
	public double getSum() {
		return getData().stream().mapToDouble(D->D.getPieValue()).sum();
	}
	
	public List<Data> getData() {
		return this.data;
		}
	@Override
	public void update() {
		
		}
	}
