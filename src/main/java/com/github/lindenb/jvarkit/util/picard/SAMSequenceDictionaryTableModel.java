/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.picard;

import javax.swing.table.AbstractTableModel;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

@SuppressWarnings("serial")
public class SAMSequenceDictionaryTableModel extends AbstractTableModel
	{
	private static final String[] COLS={"IDX","Name","Length","Species","Ideogram"};
	private SAMSequenceDictionary dict=null;
	private double longest=1L;
	public SAMSequenceDictionaryTableModel(SAMSequenceDictionary dict)
		{
		this.dict=dict;
		if(dict!=null)
			{
			for(SAMSequenceRecord ssr:dict.getSequences())
				{
				longest=Math.max(ssr.getSequenceLength(),longest);
				}
			}
		}
	
	@Override
	public int getColumnCount() {
		return COLS.length;
		}
	
	@Override
	public String getColumnName(int column) {
		return COLS[column];
		}

	@Override
	public int getRowCount()
		{
		return dict==null?0:dict.getSequences().size();
		}

	@Override
	public Object getValueAt(int row, int col)
		{
		if(dict==null) return null;
		SAMSequenceRecord ssr=dict.getSequence(row);
		if(ssr==null) return null;
		switch(col)
			{
			case 0: return ssr.getSequenceIndex();
			case 1: return ssr.getSequenceName();
			case 2: return ssr.getSequenceLength();
			case 3: return ssr.getSpecies();
			case 4: 
				{
				int n=(int)((ssr.getSequenceLength()/this.longest)*49.0);
				StringBuilder b=new StringBuilder(n+1);
				for(int i=0;i<= n;++i)
					{
					b.append("\u25A0");
					}
				return b.toString();
				}
			}
		return null;
		}
	
	@Override
	public Class<?> getColumnClass(int columnIndex)
		{
		switch(columnIndex)
			{
			case 0: return Integer.class;
			case 1: return String.class;
			case 2: return Integer.class;
			case 3: return String.class;
			case 4: return String.class;
			}
		return Object.class;
		}
	
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
		}

}
