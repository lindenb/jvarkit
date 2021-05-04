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
		case 6: return vc.isFiltered()?String.join(";",vc.getFilters()):null;
		default: throw new IllegalStateException();
		}
	}
}
