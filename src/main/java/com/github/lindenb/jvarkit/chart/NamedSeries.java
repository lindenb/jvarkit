package com.github.lindenb.jvarkit.chart;

import java.util.Collections;
import java.util.List;


public class NamedSeries extends AbstractSeries<Double> {
	public NamedSeries(final String name, List<Double> values) {
		super(name,values);
		}
	public NamedSeries(final String name,double value) {
		super(name,Collections.singletonList(value));
		}
	public double sum() {
		return stream().mapToDouble(Double::doubleValue).sum();
		}
	}
