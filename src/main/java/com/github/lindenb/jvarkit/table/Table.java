package com.github.lindenb.jvarkit.table;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public interface Table extends Iterable<List<Object>> {
public String getTitle();
public void setTitle(final String title);

public int getColumnCount();
public List<? extends Column> getColumns();
public List<List<Object>> getRows();
public default List<Object> getRow(int index) {
	return getRows().get(index);
	}
public default Column getColumn(int index) {
	return getColumns().get(index);
	}

@Override
public default Iterator<List<Object>> iterator() {
	return getRows().iterator();
	}

public default int getRowCount() {
	return getRows().size();
	}

public default boolean isEmpty() {
	return getRows().isEmpty();
	}
public void insertNewColumn(int index,final String colName);
public default void appendNewColumn(final String colName) {
	insertNewColumn(getColumnCount(), colName);
	}
public default void addRow(final Object...cells) {
	if(cells.length==1 && cells[0] instanceof List) {
		this.addRow((List<Object>)(cells[0]));
		}
	else
		{
		this.addRow(Arrays.asList(cells));
		}
	}
public void addRow(final List<Object> cells);

public default Object at(int y,int x) {
	return getRows().get(y).get(x);
}
public default void setAt(int y,int x,final Object o) {
	getRow(y).set(x, o);
	}
public void removeEmptyColumns();
}
