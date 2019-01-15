package com.github.lindenb.jvarkit.table;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class TableFactory {

private static final Function<Object, String> DEFAULT_TOSTRING = O->O==null?"":String.valueOf(O).toString();

public Table createEmptyTable() {
	return createTable(Collections.emptyList());
	}

public Table createTable(String...labels) {
	return createTable(Arrays.asList(labels));
	}

public Table createTable(final List<String> labels) {
	final TableImpl t = new TableImpl();
	t.columns.addAll(
			labels.stream().map(S->new ColumnImpl(t,S)).
			collect(Collectors.toList())
			);
	return t;
	}

private class TableImpl 
	implements Table
	{
	private String title="";
	private final List<List<Object>> rows = new ArrayList<>();
	private final List<ColumnImpl> columns = new ArrayList<>();
	private final Function<Object, String> funToString;
	
	TableImpl() {
		this.funToString = DEFAULT_TOSTRING;
		}
	
	@Override
	public void setTitle(String title) {
		this.title = title;
		}
	@Override
	public String getTitle() {
		return this.title;
		}
	@Override
	public int getColumnCount() {
		return this.columns.size();
		}
	@Override
	public List<? extends Column> getColumns() {
		return this.columns;
		}	
	@Override
	public List<List<Object>> getRows() {
		return this.rows;
		}
	
	public void insertNewColumn(int index,final String colName) {
		final ColumnImpl c = new ColumnImpl(this,colName);
		this.columns.add(index,c);
		for(int y=0;y< this.rows.size();y++) {
			this.rows.get(y).add(index, null);
			}
		}
	
	@Override
	public void addRow(final List<Object> cells) {
		if(cells.size()!=this.columns.size()) {
			throw new IllegalArgumentException("row.size:"+cells.size()+" but columns.size:"+this.columns.size());
		}
		final List<Object> r = new ArrayList<>(cells);
		this.rows.add(r);
		for(int x=0;x<this.columns.size();++x)
			{
			this.columns.get(x).maxLen=-1;
			}
		}
	@Override
	public void removeEmptyColumns()
		{
		int i=0;
		while(i < this.columns.size())
			{
			final ColumnImpl col = this.columns.get(i);
			boolean empty=true;
			for(final List<Object> row:this.rows)
				{
				if(i>=row.size()) continue;
				final String o= col.toString(row.get(i));
				if(StringUtils.isBlank(o)) continue;
				empty=false;
				break;
				}
			if(empty)
				{
				this.columns.remove(i);
				for(final List<Object> row:this.rows)
					{
					if(i>=row.size()) continue;
					row.remove(i);
					}
				}
			else
				{
				++i;
				}
			}
		}
	}

private static class ColumnImpl
	implements Column
	{
	private final TableImpl owner;
	private String name="";
	private int maxLen = 0 ;
	private Function<Object, String> funToString;
	ColumnImpl(final TableImpl owner,final String name) {
		this.owner = owner;
		this.funToString  = owner.funToString;
		this.name = StringUtils.isBlank(name)?"":name;
		this.maxLen = this.name.length();
		}
	
	@Override
	public void setStringConverter(final Function<Object, String> fun) {
		this.funToString = fun;
		}
	
	@Override
	public String toString(final Object o) {
		final String s;
		if(this.funToString!=null)
			{
			s = this.funToString.apply(o);
			}
		else if(this.owner.funToString!=null)
			{
			s = this.owner.funToString.apply(o);
			}
		else
			{
			s = DEFAULT_TOSTRING.apply(o);
			}
		return s==null?"":s;
		}
	
	public String getName() {
		return this.name;
		}
	@Override
	public String toString() {
		return getName();
		}
	}


}
