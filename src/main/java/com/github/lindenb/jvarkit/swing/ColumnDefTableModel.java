/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.swing;

import java.util.Collections;
import java.util.List;
import java.util.Vector;


public class ColumnDefTableModel<DATATYPE> extends AbstractGenericTableModel<DATATYPE> {
	private static final long serialVersionUID = 1L;
	private final List<ColumnDef<DATATYPE>> columns = new Vector<>();
	
	public ColumnDefTableModel(final List<ColumnDef<DATATYPE>> columns) {
		this(columns,Collections.emptyList());
		}
	public ColumnDefTableModel(final List<ColumnDef<DATATYPE>> columns,final List<DATATYPE> rows) {
		super(rows);
		this.columns.addAll(columns);
		}
	public void addColumn(final ColumnDef<DATATYPE> def) {
		this.columns.add(def);
		fireTableStructureChanged();
		}	
	@Override
	public Object getValueOf(DATATYPE o, int columnIndex) {
		return this.columns.get(columnIndex).getExtractor().apply(o);
		}
	
	public Class<?> getColumnClass(int column) {
		return this.columns.get(column).getColumnClass();
		}
	public String getColumnName(int column) {
		return this.columns.get(column).getColumnName();
		}
	@Override
	public final int getColumnCount() {
		return columns.size();
		}
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
		}
	}
