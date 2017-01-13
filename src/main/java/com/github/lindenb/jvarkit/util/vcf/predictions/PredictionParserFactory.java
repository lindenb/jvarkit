package com.github.lindenb.jvarkit.util.vcf.predictions;

import htsjdk.variant.vcf.VCFHeader;

public class PredictionParserFactory {
	private String tag="ANN";
	public void setTag(String tag) {
		this.tag = tag;
	}
	public String getTag() {
		return tag;
	}
	
	public AnnPredictionParser buildAnnPredictionParser(final VCFHeader header) {
		return new AnnPredictionParser(header, getTag());
	}
}
