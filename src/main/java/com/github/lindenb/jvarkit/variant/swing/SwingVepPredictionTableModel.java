package com.github.lindenb.jvarkit.variant.swing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingVepPredictionTableModel extends AbstractGenericTable<VepPrediction>
	{
	private final List<String> COLS;

	private final VepPredictionParser vepPredictionParser;
	public SwingVepPredictionTableModel(final VCFHeader header) {
		this.vepPredictionParser = new VepPredictionParserFactory(header).get();
		this.COLS=new ArrayList<>(this.vepPredictionParser.getCategories());
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
	
	public void setVariant(final VariantContext vc) {
		final List<VepPrediction> L;
		if(vc!=null && vc.hasAttribute(this.vepPredictionParser.getTag())) {
			L = this.vepPredictionParser.getPredictions(vc);
			}
		else
			{
			L = Collections.emptyList();
			}
		this.setRows(L);
		}
		
	@Override
	public Object getValueOf(final VepPrediction P, int columnIndex)
		{
		return P.get(this.getColumnName(columnIndex));
		}
	}