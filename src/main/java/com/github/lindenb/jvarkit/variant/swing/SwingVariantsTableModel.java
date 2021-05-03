package com.github.lindenb.jvarkit.variant.swing;

import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("serial")
public class SwingVariantsTableModel extends AbstractGenericTable<VariantContext>{
public SwingVariantsTableModel() {
	super();
	}

@Override
public int getColumnCount() {
	return 7;
	}

@Override
public Class<?> getColumnClass(int column) {
		switch(column) {
		case 0: return String.class;
		case 1: return Integer.class;
		case 2: return String.class;
		case 3: return String.class;
		case 4: return String.class;
		case 5: return Double.class;
		case 6: return String.class;
		default: throw new IllegalStateException();
		}
	}

@Override
public String getColumnName(int column) {
	switch(column) {
		case 0: return "CHROM";
		case 1: return "POS";
		case 2: return "ID";
		case 3: return "REF";
		case 4: return "ALT";
		case 5: return "QUAL";
		case 6: return "FILTER";
		default: throw new IllegalStateException();
		}
	}

@Override
public Object getValueOf(final VariantContext vc, int column) {
		switch(column) {
		case 0: return vc.getContig();
		case 1: return vc.getStart();
		case 2: return vc.hasID()?vc.getID():null;
		case 3: return vc.getReference().getDisplayString();
		case 4: return vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","));
		case 5: return vc.hasLog10PError()?vc.getPhredScaledQual():null;
		case 6: return vc.isFiltered()?String.join(",",vc.getFilters()):null;
		default: throw new IllegalStateException();
		}
	}
}
