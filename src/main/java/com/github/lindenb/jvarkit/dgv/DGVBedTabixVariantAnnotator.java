/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.dgv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/*** tabix annotator using DGV SV as BED */
public class DGVBedTabixVariantAnnotator extends AbstractTabixVariantAnnotator {
	public static final String OPT_DESC="DGV SV Variant file as Tabix indexed file from http://dgv.tcag.ca/dgv/app/downloads";
	public static final String FRACTION_KEY ="fraction";
	private final FileHeader bedHeader;
	private double fraction;
	private VCFInfoHeaderLine hdrHasDGV;
	private VCFInfoHeaderLine hdrDGVSV;
	private VCFInfoHeaderLine hdrDGVAF;
	
	public DGVBedTabixVariantAnnotator(final String bedUri) throws IOException {
		this(bedUri,AttributeMap.empty());
		}
	
	/**
	 * 
	 * @param bedUri
	 * @param meta 'fraction': common overlap fraction
	 * @throws IOException
	 */
	public DGVBedTabixVariantAnnotator(final String bedUri,  AttributeMap meta) throws IOException {
		super(bedUri);
		if(super.tabixReader!=null) {
			if(meta==null) meta = AttributeMap.empty();
			final String line = super.tabixReader.readLine();
			if(line==null) throw new IOException("Cannot read header from "+bedUri);
			this.bedHeader = new FileHeader(CharSplitter.TAB.split(line));
			this.fraction = meta.getDoubleAttribute(FRACTION_KEY).orElse(0.9);
			}
		else
			{
			this.bedHeader = null;
			this.fraction = 0.0;
			}
		}
	
	/** set fraction [0-1] of required multual overlap */
	public void setFraction(double fraction) {
		this.fraction = fraction;
		}
	
	public double getFraction() {
		return fraction;
		}
	
	@Override
	public void fillHeader(VCFHeader header) {
		if(super.tabixReader==null) return;
		final String cmdHelp = "DGV "+super.getUri()+" with shared overlap of "+this.getFraction();
		this.hdrHasDGV = new VCFInfoHeaderLine("HAS_DGV_SV",1, VCFHeaderLineType.Flag, "SV was found in "+cmdHelp);
		this.hdrDGVSV = new VCFInfoHeaderLine("DGV_SV",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, 
				"SV in "+ cmdHelp +" Format:ID|SVTYPE|Frequency");
		this.hdrDGVAF = new VCFInfoHeaderLine("DGV_FREQ_MAX",1, VCFHeaderLineType.Float, "Max frequency found in " + cmdHelp);

		header.addMetaDataLine(this.hdrHasDGV);
		header.addMetaDataLine(this.hdrDGVSV);
		header.addMetaDataLine(this.hdrDGVAF);
		}
	
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx)  throws IOException {
		if(super.tabixReader==null || !hasContig(ctx)) return Collections.singletonList(ctx);
		final Set<String> dgv_variants = new HashSet<>();
		TabixReader.Iterator r = super.tabixReader.query(contig(ctx),ctx.getStart()-1, ctx.getEnd());
		Double popmax_af = null;
		for(;;) {
			String line = r.next();
			if(line==null) break;
			final Map<String, String> rec = this.bedHeader.toMap(CharSplitter.TAB.split(line));
			final int svStart = Integer.parseInt(rec.get("start"))+1;
			final int svEnd = Integer.parseInt(rec.get("end"));
			if(!this.overlaps(ctx, svStart, svEnd,getFraction())) continue;
			
			Double rec_af = null;
			final OptionalInt sample_size = StringUtils.parseInt(rec.get("samplesize"));
			if(sample_size.isPresent() && sample_size.getAsInt()>0) {
				final OptionalInt observedgains = StringUtils.parseInt(rec.get("observedgains"));
				final OptionalInt observedlosses = StringUtils.parseInt(rec.get("observedlosses"));
				rec_af = (observedgains.orElse(0)+observedlosses.orElse(0))/(double)sample_size.getAsInt();
				}
			if(popmax_af==null || popmax_af.doubleValue()< rec_af.doubleValue()) {
				popmax_af = rec_af.doubleValue();
				}
			String type = rec.get("variantsubtype");
			if(StringUtils.isBlank(type)) type= rec.get("varianttype");
			if(StringUtils.isBlank(type)) type= "";
			type = type.replace(' ', '_');//eg 'line deletion'
			dgv_variants.add(String.join("|", rec.get("variantaccession"),type,(rec_af!=null?rec_af.toString():".")));
			}

		if(!dgv_variants.isEmpty()) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.attribute(this.hdrHasDGV.getID(), true);
			vcb.attribute(this.hdrDGVSV.getID(), new ArrayList<>(dgv_variants));
			if(popmax_af!=null) vcb.attribute(this.hdrDGVAF.getID(), popmax_af);
			return Collections.singletonList(vcb.make());
			}
		else
			{
			return Collections.singletonList(ctx);
			}
		}
	
}
