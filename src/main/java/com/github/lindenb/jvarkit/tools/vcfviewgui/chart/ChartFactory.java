package com.github.lindenb.jvarkit.tools.vcfviewgui.chart;

import java.util.Collection;

import com.github.lindenb.jvarkit.tools.vcfviewgui.PedFile;

import javafx.scene.chart.Chart;

public interface ChartFactory<HEADER,T> {
public String getName();
public void setHeader(final HEADER header);
public HEADER getHeader();
public PedFile getPedigree();
public void setPedigree(final PedFile pedigree);
public void visit(final T o);
public default void visit(final Collection<T> L) {
	for(final T o:L) this.visit(o);
	}
public Chart build();
}
