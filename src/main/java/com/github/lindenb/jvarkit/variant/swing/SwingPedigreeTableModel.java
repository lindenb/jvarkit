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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

@SuppressWarnings("serial")
public class SwingPedigreeTableModel extends AbstractGenericTable<Sample>
	{
	private final List<String> COLS = Arrays.asList("FAM","ID","FATHER","MOTHER","SEX","STATUS");

	public SwingPedigreeTableModel(final Pedigree ped) {
		super(new ArrayList<>(ped.getSamples()));
	}
	
	@Override
	public int getColumnCount() {
		return COLS.size();
		}
	@Override
	public String getColumnName(int column)
		{
		return COLS.get(column);
		}
	
		
	@Override
	public Object getValueOf(final Sample P, int columnIndex)
		{
		switch(columnIndex) {
			case 0: return P.getFamily();
			case 1: return P.getId();
			case 2: return P.hasFather()?P.getFather().getId():null;
			case 3: return P.hasMother()?P.getMother().getId():null;
			case 4: return P.getSex().name();
			case 5: return P.isStatusSet()?P.getStatus().name():null;
			default: throw new IllegalArgumentException();
			}
		}
	}