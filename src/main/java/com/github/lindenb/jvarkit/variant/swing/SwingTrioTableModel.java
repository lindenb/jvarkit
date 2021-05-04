package com.github.lindenb.jvarkit.variant.swing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.Trio;
import com.github.lindenb.jvarkit.tools.vcftrios.DeNovoDetector;
import com.github.lindenb.jvarkit.tools.vcftrios.DeNovoDetector.DeNovoMutation;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("serial")
public class SwingTrioTableModel extends AbstractGenericTable<DeNovoDetector.DeNovoMutation>
	{
	private final List<Trio> trios;
	private final List<String> COLS = Arrays.asList(
			"FATHER","FATHER-GT",
			"MOTHER","MOTHER-GT",
			"CHILD","CHILD-GT"
			);
	private final DeNovoDetector detector = new DeNovoDetector();
	public SwingTrioTableModel(final Pedigree ped) {
		this.trios = (ped==null?Collections.emptyList(): new ArrayList<>(ped.getTrios()));
	}
	
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return String.class;
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
	public void setVariant(final VariantContext ctx) {
		final List<DeNovoDetector.DeNovoMutation> L;
		if(ctx==null || !ctx.hasGenotypes() || this.trios.isEmpty()) {
			L = Collections.emptyList();
			}
		else {
			L = new ArrayList<>();
			for(final Trio trio:this.trios) {
				final DeNovoMutation mut = this.detector.test(ctx,trio);
				if(mut!=null) L.add(mut);
				}
			}
		setRows(L);
		}
		
	@Override
	public Object getValueOf(final DeNovoMutation P, int columnIndex)
		{
		switch(columnIndex) {
			case 0: return P.hasFatherGenotype()?P.getFatherGenotype().getSampleName():null;
			case 1: return P.hasFatherGenotype()?P.getFatherGenotype().getGenotypeString():null;
			case 2: return P.hasMotherGenotype()?P.getMotherGenotype().getSampleName():null;
			case 3: return P.hasMotherGenotype()?P.getMotherGenotype().getGenotypeString():null;
			case 4: return P.hasChildGenotype()?P.getChildGenotype().getSampleName():null;
			case 5: return P.hasChildGenotype()?P.getChildGenotype().getGenotypeString():null;
			default: throw new IllegalArgumentException();
			}
		}
	}