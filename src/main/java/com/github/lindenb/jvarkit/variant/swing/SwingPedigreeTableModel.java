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