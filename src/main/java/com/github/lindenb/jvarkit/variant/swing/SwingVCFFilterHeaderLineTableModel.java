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
