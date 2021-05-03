package com.github.lindenb.jvarkit.variant.swing;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Vector;

import javax.swing.table.AbstractTableModel;

import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("serial")
public class SwingVCFInfoTableModel extends AbstractTableModel {
private final Vector<InfoLine> rows = new Vector<>();
private class InfoLine {
	String key;
	Integer idx;
	Object value=null;
	}

public SwingVCFInfoTableModel() {
	super();
	}

public void setVariant(final VariantContext ctx) {
	setAttributes(ctx==null?null:ctx.getAttributes());
	}

public void setAttributes(final Map<String,Object> map) {
	this.rows.clear();
	if(map!=null && !map.isEmpty()) {
		for(final String key: new TreeSet<>(map.keySet()))
			{
			Object v= map.get(key);
			final List<?> L;
			if(v instanceof List)
				{
				L=(List<?>)v;
				}
			else if(v.getClass().isArray())
				{
				Object a[]=(Object[])v;
				L=Arrays.asList(a);
				}
			else
				{
				L=Collections.singletonList(v);
				}
			for(int x=0;x< L.size();++x)
				{
				final InfoLine F = new InfoLine();
				F.key = key;
				F.idx = (L.size()==1?null:x+1);
				F.value = L.get(x);
				rows.add(F);
				}
			}
		}
	
	fireTableDataChanged();
	}


@Override
public int getColumnCount() {
	return 3;
	}

@Override
public int getRowCount() {
		return rows.size();
	}

@Override
public Class<?> getColumnClass(int column) {
	switch(column) {
		case 0: return String.class;
		case 1: return Integer.class;
		case 2: return Object.class;
		default: throw new IllegalArgumentException();
		}
	}
@Override
public String getColumnName(int column) {
	switch(column) {
	case 0: return "ID";
	case 1: return "Index";
	case 2: return "Value";
	default: throw new IllegalArgumentException();
	}
}
@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
	}

@Override
public Object getValueAt(int rowIndex, int column) {
	final InfoLine L = this.rows.get(rowIndex);
	switch(column) {
	case 0: return L.key;
	case 1: return L.idx;
	case 2: return L.value;
	default: throw new IllegalArgumentException();
	}
	}
}
