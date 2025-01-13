/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.gnomad;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;

import htsjdk.samtools.util.CoordMath;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/*** tabix annotator using gnomad SV as BED */
public class GnomadSVBedTabixVariantAnnotator extends AbstractTabixVariantAnnotator {
	public static final String OPT_DESC="Gnomad SV as BED file. Compressed with bgzip and indexed with tabix.";
	public static final String FRACTION_KEY ="fraction";

	private final FileHeader bedHeader;
	private final double fraction;
	private VCFInfoHeaderLine hdrHasGnomad;
	private VCFInfoHeaderLine hdrGnomadSV;
	private VCFInfoHeaderLine hdrGnomadMaxAF;
	private final String af_column;
	/**
	 * 
	 * @param bedUri
	 * @param meta map with 'af_column': AF column to watch, 'fraction': common overlap fraction
	 * @throws IOException
	 */
	public GnomadSVBedTabixVariantAnnotator(final String bedUri,  AttributeMap meta) throws IOException {
		super(bedUri);
		if(super.tabixReader!=null) {
			if(meta==null) meta = AttributeMap.empty();
			final String line = super.tabixReader.readLine();
			if(line==null) throw new IOException("Cannot read header from "+bedUri);
			this.bedHeader = new FileHeader(CharSplitter.TAB.split(line));
			this.fraction = meta.getDoubleAttribute(FRACTION_KEY).orElse(0.9);
			this.af_column = meta.getAttribute("af_column").orElse("AF");
			this.bedHeader.getColumnIndex(this.af_column);
			}
		else
			{
			this.bedHeader = null;
			this.af_column = null;
			this.fraction = 0.0;
			}
		}
	
	@Override
	public void fillHeader(VCFHeader header) {
		if(super.tabixReader==null) return;
		final String cmdHelp = "gnomadSV "+super.getUri()+" with shared overlap of "+this.fraction+".";
		this.hdrHasGnomad = new VCFInfoHeaderLine("HAS_GNOMAD_SV",1, VCFHeaderLineType.Flag, "SV was found in "+cmdHelp);
		this.hdrGnomadSV = new VCFInfoHeaderLine("GNOMAD_SV",
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String, 
				"SV in "+ cmdHelp +" Format:ID|SVTYPE|"+af_column);
		this.hdrGnomadMaxAF = new VCFInfoHeaderLine("GNOMAD_SV_"+this.af_column+"_MAX",1,
				VCFHeaderLineType.Float,
				"Max "+this.af_column+" frequency found in " + cmdHelp
				);

		header.addMetaDataLine(this.hdrHasGnomad);
		header.addMetaDataLine(this.hdrGnomadSV);
		header.addMetaDataLine(this.hdrGnomadMaxAF);
		}
	@Override
	public List<VariantContext> annotate(final VariantContext ctx)  throws IOException {
		if(super.tabixReader==null || !hasContig(ctx)) return Collections.singletonList(ctx);
		final Set<String> gnomad_variants = new HashSet<>();
		TabixReader.Iterator r = super.tabixReader.query(contig(ctx),ctx.getStart()-1, ctx.getEnd());
		Double popmax_af = null;
		for(;;) {
			final String line = r.next();
			if(line==null) break;
			final Map<String, String> rec = this.bedHeader.toMap(CharSplitter.TAB.split(line));
			final int svStart = Integer.parseInt(rec.get("start"))+1;
			final int svEnd = Integer.parseInt(rec.get("end"));
			if(!CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), svStart, svEnd)) continue;
			
			final int x0 = Math.max(ctx.getStart(),svStart);
			final int x1 = Math.min(ctx.getEnd(),svEnd);
			final double shared_len = CoordMath.getLength(x0, x1);
			if(shared_len/ctx.getLengthOnReference() < this.fraction) continue; 
			if(shared_len/CoordMath.getLength(svStart, svEnd) < this.fraction) continue;
			final OptionalDouble rec_popmax_af = Arrays.stream(CharSplitter.COMMA.split(rec.get(this.af_column))).
					filter(S->!S.equals("NA")).
					mapToDouble(Double::parseDouble).
					max()
					;
			gnomad_variants.add(String.join("|", rec.get("name"),rec.get("svtype"),(rec_popmax_af.isPresent()?String.valueOf(rec_popmax_af.getAsDouble()):".")));
			if(rec_popmax_af.isPresent()) {
				if(popmax_af==null || popmax_af.doubleValue()< rec_popmax_af.getAsDouble()) {
					popmax_af = rec_popmax_af.getAsDouble();
					}
				}
			
			}
		if(!gnomad_variants.isEmpty()) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.attribute(this.hdrHasGnomad.getID(), true);
			vcb.attribute(this.hdrGnomadSV.getID(), new ArrayList<>(gnomad_variants));
			if(popmax_af!=null) vcb.attribute(this.hdrGnomadMaxAF.getID(), popmax_af);
			return Collections.singletonList(vcb.make());
			}
		else
			{
			return Collections.singletonList(ctx);
			}
		}
	
}
