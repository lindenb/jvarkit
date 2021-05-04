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

import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

@SuppressWarnings("serial")
public class SwingAnnPredictionTableModel extends AbstractGenericTable<AnnPrediction>
	{
	private final String[] COLS = new String[]{
			"SO","Allele","Impact","GeneName","GeneId","FeatureType","FeatureId",
			"BioType","HGVsc","HGVsp","Rank","cDNA-pos","CDS-pos","AA-pos","Distance","Msg"};

	private final AnnPredictionParser annPredictionParser;
	public SwingAnnPredictionTableModel(final VCFHeader header) {
		this.annPredictionParser = new AnnPredictionParserFactory().header(header).get();
	}
	
	@Override
	public int getColumnCount() {
		return COLS.length;
		}
	@Override
	public String getColumnName(int column)
		{
		return COLS[column];
		}
	
	public void setVariant(final VariantContext vc) {
		final List<AnnPrediction> L;
		if(vc!=null && vc.hasAttribute(this.annPredictionParser.getTag())) {
			L = this.annPredictionParser.getPredictions(vc);
			}
		else
			{
			L = Collections.emptyList();
			}
		this.setRows(L);
		}
		
	@Override
	public Object getValueOf(AnnPrediction P, int columnIndex)
		{
		switch(columnIndex) {
			case 0: return P.getSOTermsString();
			case 1: return P.getAllele();
			case 2: return P.getPutativeImpact();
			case 3: return P.getGeneName();
			case 4: return P.getGeneId();
			case 5: return P.getFeatureType();
			case 6: return P.getFeatureId();
			case 7: return P.getTranscriptBioType();
			case 8: return P.getHGVSc();
			case 9: return P.getHGVSp();
			case 10: return P.getRank();
			case 11: return P.getCDNAPos();
			case 12: return P.getCDSPos();
			case 13: return P.getAAPos();
			case 14: return P.getDistance();
			case 15: return P.getMessages();
		    default: throw new IllegalArgumentException();
			}			
		}
	}