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
package com.github.lindenb.jvarkit.gff3;

import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.tribble.gff.Gff3Feature;

@SuppressWarnings("serial")
public class SwingGff3TableModel extends AbstractGenericTable<Gff3Feature>{
private final Vector<ColumnInfo> columns = new Vector<>();
private class ColumnInfo {
	String name;
	Class<?> clazz;
	Function<Gff3Feature, Object> extractor;
	}

public SwingGff3TableModel() {
	}


@Override
synchronized public void  setRows(final List<Gff3Feature> features) {
	if(features==null || features.isEmpty()) {
		setRows(Collections.emptyList());
		}
    this.columns.clear();
    
    final Vector<ColumnInfo> columns = new Vector<>();
    
    ColumnInfo ci = new ColumnInfo();
    ci.clazz = String.class;
    ci.name= "Contig";
    ci.extractor = GT->GT.getContig();
    columns.add(ci);
    
    ci = new ColumnInfo();
    ci.clazz = String.class;
    ci.name= "Source";
    ci.extractor = GT->GT.getSource();
    columns.add(ci);
    
    ci = new ColumnInfo();
    ci.clazz = String.class;
    ci.name= "Type";
    ci.extractor = GT->GT.getType();
    columns.add(ci);

    ci = new ColumnInfo();
    ci.clazz = Integer.class;
    ci.name= "Start";
    ci.extractor = GT->GT.getStart();
    columns.add(ci);
    
    ci = new ColumnInfo();
    ci.clazz = Integer.class;
    ci.name= "End";
    ci.extractor = GT->GT.getEnd();
    columns.add(ci);
    
    ci = new ColumnInfo();
    ci.clazz = Integer.class;
    ci.name= "Phase";
    ci.extractor = GT->GT.getPhase();
    columns.add(ci);

    
    ci = new ColumnInfo();
    ci.clazz = String.class;
    ci.name= "Strand";
    ci.extractor = GT->GT.getStrand()==null?null:GT.getStrand().name();
    columns.add(ci);
    
    ci = new ColumnInfo();
    ci.clazz = Double.class;
    ci.name= "Score";
    ci.extractor = GT->GT.getScore();
    columns.add(ci);
    
    for(final String key:features.stream().
    	flatMap(F->F.getAttributes().keySet().stream()).collect(Collectors.toSet()))
    	{
	    ci = new ColumnInfo();
	    ci.clazz = String.class;
	    ci.name= key;
	    ci.extractor = GT->GT.getAttribute(key).stream().collect(Collectors.joining(" "));
	    columns.add(ci);
    	}
    
    this.columns.addAll(columns);
	
    
    super.rows.clear();
    super.rows.addAll(features);
	fireTableStructureChanged();
	}

@Override
public int getColumnCount() {
	return this.columns.size();
	}

@Override
public Class<?> getColumnClass(int column) {
	return this.columns.get(column).clazz;
	}
@Override
public String getColumnName(int column) {
	return this.columns.get(column).name;
	}

@Override
public Object getValueOf(final Gff3Feature F, int column) {
	return this.columns.get(column).extractor.apply(F); 
	}
}
