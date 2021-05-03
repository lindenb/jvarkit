package com.github.lindenb.jvarkit.variant.swing;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingVCFFormatHeaderLineTableModel extends AbstractGenericTable<VCFFormatHeaderLine>{
public SwingVCFFormatHeaderLineTableModel(final List<VCFFormatHeaderLine> list) {
	super(list);
	}

public SwingVCFFormatHeaderLineTableModel(final VCFHeader header) {
	this(header.getFormatHeaderLines().stream().
			sorted((A,B)->A.getID().compareTo(B.getID())).
			collect(Collectors.toList()));
	}

@Override
public int getColumnCount() {
	return 5;
	}

@Override
public Class<?> getColumnClass(int column) {
	switch(column) {
		case 2: return Integer.class;
		default:return String.class;
		}
	}
@Override
public String getColumnName(int column) {
	switch(column) {
		case 0: return "ID";
		case 1: return "Type";
		case 2: return "Count";
		case 3: return "CountType";
		case 4: return "Description";
		default: throw new IllegalStateException();
		}
	}

@Override
public Object getValueOf(final VCFFormatHeaderLine F, int columnIndex) {
	switch(columnIndex) {
		case 0: return F.getID();
		case 1 : return F.getType()==null?null:F.getType().name();
		case 2: return F.isFixedCount()?F.getCount():null;
		case 3: return F.getCountType().name();
		case 4: return F.getDescription();
		default: throw new IllegalStateException();
		}
	}
}
