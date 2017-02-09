package com.github.lindenb.jvarkit.tools.vcfviewgui.chart;

import java.util.Collection;

import javafx.scene.chart.Chart;

public interface ChartFactory<T> {
public String getName();
public void visit(final T o);
public default void visit(final Collection<T> L) {
	for(final T o:L) this.visit(o);
	}
public Chart build();
}
