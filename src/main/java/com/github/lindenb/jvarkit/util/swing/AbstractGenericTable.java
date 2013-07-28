package com.github.lindenb.jvarkit.util.swing;

import java.util.List;
import java.util.Vector;
import javax.swing.table.AbstractTableModel;


@SuppressWarnings("serial")
public abstract class AbstractGenericTable<T>
	extends AbstractTableModel
	{
	protected List<T> rows=null;
	
	
	public AbstractGenericTable()
		{
		this.rows=new Vector<T>();
		}
	
	public AbstractGenericTable(List<T> rows)
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
	
	public synchronized void setRows(List<T> L)
		{
		this.rows=new Vector<T>(L);
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
