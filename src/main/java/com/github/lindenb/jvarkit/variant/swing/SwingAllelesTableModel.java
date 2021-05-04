/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.variant.swing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.table.AbstractTableModel;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("serial")
public class SwingAllelesTableModel extends AbstractTableModel {
private final List<Allele> alleles = new ArrayList<>();

public void setVariant(final VariantContext ctx) {
	setRows(ctx==null?Collections.emptyList():ctx.getAlleles());
	}

public void setRows(final List<Allele> alleles) {
	this.alleles.clear();
	if(alleles!=null) this.alleles.addAll(alleles);
	fireTableDataChanged();
	}


@Override
public int getColumnCount() {
	return 5;
	}

@Override
public int getRowCount() {
		return this.alleles.size();
	}

@Override
public Class<?> getColumnClass(int column) {
	switch(column) {
		case 0: return Integer.class;
		case 1: return String.class;
		case 2: return Boolean.class;
		case 3: return Boolean.class;
		case 4: return Integer.class;
		default: throw new IllegalArgumentException();
		}
	}
@Override
public String getColumnName(int column) {
	switch(column) {
	case 0: return "Index";
	case 1: return "Bases";
	case 2: return "REF";
	case 3: return "Symbolic";
	case 4: return "Length";
	default: throw new IllegalArgumentException();
	}
}
@Override
public boolean isCellEditable(int rowIndex, int columnIndex) {
	return false;
	}

@Override
public Object getValueAt(int rowIndex, int column) {
	final Allele a = this.alleles.get(rowIndex);
	switch(column) {
	case 0: return rowIndex;
	case 1: return a.getDisplayString();
	case 2: return a.isReference();
	case 3: return a.isSymbolic();
	case 4: return a.length();
	default: throw new IllegalArgumentException();
	}
	}
}
