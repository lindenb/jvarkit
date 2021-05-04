package com.github.lindenb.jvarkit.variant.swing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser.BcfToolsPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingBcsqPredictionTableModel extends AbstractGenericTable<BcfToolsPrediction>
	{
	private final List<String> COLS;

	private final BcfToolsPredictionParser bcfPredictionParser;
	public SwingBcsqPredictionTableModel(final VCFHeader header) {
		this.bcfPredictionParser = new BcfToolsPredictionParserFactory(header).get();
		this.COLS=new ArrayList<>(this.bcfPredictionParser.getCategories());
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
		final List<BcfToolsPrediction> L;
		if(vc!=null && vc.hasAttribute(this.bcfPredictionParser.getTag())) {
			L = this.bcfPredictionParser.getPredictions(vc);
			}
		else
			{
			L = Collections.emptyList();
			}
		this.setRows(L);
		}
		
	@Override
	public Object getValueOf(final BcfToolsPrediction P, int columnIndex)
		{
		return P.getByCol(this.getColumnName(columnIndex));
		}
	}