/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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

import java.util.List;
import java.util.Vector;

import javax.swing.table.AbstractTableModel;


@SuppressWarnings("serial")
public abstract class AbstractGenericTableModel<T>
	extends AbstractTableModel
	{
	protected final List<T> rows;
	
	public AbstractGenericTableModel()
		{
		this.rows=new Vector<T>();
		}
	
	public AbstractGenericTableModel(final List<T> rows)
		{
		this.rows=new Vector<T>(rows);
		}
	
	public List<T> getRows()
		{
		return this.rows;
		}
	
	public void clear()
		{
		this.rows.clear();
		fireTableDataChanged();
		}
	
	public synchronized void setRows(final List<T> L)
		{
		this.rows.clear();
		this.rows.addAll(L);
		fireTableDataChanged();
		}
	
	@Override
	public abstract int getColumnCount();

	@Override
	public final int getRowCount()
		{
		return getRows().size();
		}
	
	public abstract Object getValueOf(T o, int columnIndex);
	
	@Override
	public final Object getValueAt(int rowIndex, int columnIndex)
		{
		T o=getElementAt(rowIndex);
		if(o==null) return null;
		return getValueOf(o, columnIndex);
		}
	
	public T getElementAt(int rowIndex)
		{
		return  getRows().get(rowIndex);
		}
	
	
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex)
		{
		return false;
		}
	}
