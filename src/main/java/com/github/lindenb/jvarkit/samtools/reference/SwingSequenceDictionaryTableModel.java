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
package com.github.lindenb.jvarkit.samtools.reference;

import java.util.List;
import java.util.Vector;
import java.util.stream.Collectors;

import javax.swing.table.AbstractTableModel;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

@SuppressWarnings("serial")
public class SwingSequenceDictionaryTableModel extends AbstractTableModel {
	final SAMSequenceDictionary dict;
	final List<String> keys;
	
	public SwingSequenceDictionaryTableModel(final SAMSequenceDictionary dict) {
		this.dict = dict;
		this.keys=new Vector<>(dict.getSequences().stream().
				flatMap(SSR->SSR.getAttributes().stream()).
				map(KV->KV.getKey()).
				collect(Collectors.toSet()));
		}
	public SAMSequenceDictionary getSequenceDictionary() { return dict;}
	
	@Override
	public String getColumnName(int column)
		{
		switch (column)
			{
			case 0: return "SEQ";
			case 1: return "LENGTH";
			default: return this.keys.get(column-2);
			}
		}
	
	
	@Override
	public int getColumnCount()
		{
		return 2+this.keys.size();
		}
	@Override
	public int getRowCount()
		{
		return getSequenceDictionary().size();
		}
	@Override
	public Class<?> getColumnClass(int column)
		{
		switch (column)
			{
			case 0: return String.class;
			case 1: return Integer.class;
			default: return super.getColumnClass(column);
			}
		}
	@Override
	public Object getValueAt(int rowIndex, int column)
		{
		final SAMSequenceRecord ssr = getSequenceDictionary().getSequence(rowIndex);
		switch (column)
			{
			case 0: return ssr.getSequenceName();
			case 1: return ssr.getSequenceLength();
			default: return ssr.getAttribute(this.keys.get(column-2));
			}
		}
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex)
		{
		return false;
		}
	}
