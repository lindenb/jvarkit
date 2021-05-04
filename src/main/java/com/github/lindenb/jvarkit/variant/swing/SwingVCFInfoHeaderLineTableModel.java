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

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@SuppressWarnings("serial")
public class SwingVCFInfoHeaderLineTableModel extends AbstractGenericTable<VCFInfoHeaderLine>{
public SwingVCFInfoHeaderLineTableModel(final List<VCFInfoHeaderLine> list) {
	super(list);
	}

public SwingVCFInfoHeaderLineTableModel(final VCFHeader header) {
	this(header.getInfoHeaderLines().stream().
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
public Object getValueOf(final VCFInfoHeaderLine F, int columnIndex) {
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
