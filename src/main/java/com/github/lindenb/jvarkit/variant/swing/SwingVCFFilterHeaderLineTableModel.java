package com.github.lindenb.jvarkit.variant.swing;

import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingVCFFilterHeaderLineTableModel extends AbstractGenericTable<VCFFilterHeaderLine>{
public SwingVCFFilterHeaderLineTableModel(final List<VCFFilterHeaderLine> list) {
	super(list);
	}

public SwingVCFFilterHeaderLineTableModel(final VCFHeader header) {
	this(header.getFilterLines().stream().
			sorted((A,B)->A.getID().compareTo(B.getID())).
			collect(Collectors.toList()));
	}


@Override
public int getColumnCount() {
	return 2;
	}

@Override
public Class<?> getColumnClass(int columnIndex) {
	return String.class;
	}
@Override
public String getColumnName(int column) {
	switch(column) {
		case 0: return "ID";
		case 1: return "Description";
		default: throw new IllegalStateException();
		}
	}
@Override
public Object getValueOf(final VCFFilterHeaderLine o, int columnIndex) {
	switch(columnIndex) {
		case 0: return o.getID();
		case 1: return o.getDescription();
		default: throw new IllegalStateException();
		}
	}

}
