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
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
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