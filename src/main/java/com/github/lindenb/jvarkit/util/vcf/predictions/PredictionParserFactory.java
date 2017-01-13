package com.github.lindenb.jvarkit.util.vcf.predictions;

import htsjdk.variant.vcf.VCFHeader;

public class PredictionParserFactory {
	private String tag="ANN";
	private VCFHeader vcfHeader=null;
	
	public VCFHeader getVcfHeader() {
		return vcfHeader;
	}
	
	public void setVcfHeader(final VCFHeader vcfHeader) {
		this.vcfHeader = vcfHeader;
	}
	
	public PredictionParserFactory header(final VCFHeader vcfHeader) {
		this.setVcfHeader(vcfHeader);
		return this;
	}
	
	public void setTag(String tag) {
		this.tag = tag;
	}
	public String getTag() {
		return tag;
	}
	
	public AnnPredictionParser buildAnnPredictionParser() {
		return new AnnPredictionParser(getVcfHeader(), getTag());
	}
}
