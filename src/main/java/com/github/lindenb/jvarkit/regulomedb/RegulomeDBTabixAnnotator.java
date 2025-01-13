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
package com.github.lindenb.jvarkit.regulomedb;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class RegulomeDBTabixAnnotator  extends AbstractTabixVariantAnnotator {
	public static final String OPT_DESC="RegulomeDB bed sorted, bgzipped and indexed with tabix.";
	public static final String REGEX_KEY="regex";
	public static final String EXTEND_KEY="extend";
	private final FileHeader fileHeader;
	private final Pattern acceptRegex;
	private final int extend;
	private VCFInfoHeaderLine hdrRegulomeDbMeanScore = null;
	private VCFInfoHeaderLine hdrRegulomeDbMaxScore = null;
	
	public RegulomeDBTabixAnnotator(final String bedUri,  AttributeMap meta) throws IOException {
		super(bedUri);
		if(super.tabixReader!=null) {
			if(meta==null) meta = AttributeMap.empty();
			final String line = super.tabixReader.readLine();
			if(line==null) throw new IOException("Cannot read header from "+bedUri);
			this.fileHeader = new FileHeader(CharSplitter.TAB.split(line));
			this.fileHeader.getColumnIndex("ranking");
			this.fileHeader.getColumnIndex("probability_score");
			
			if(!StringUtils.isBlank(meta.getAttribute(REGEX_KEY, "")))
				{
				this.acceptRegex=Pattern.compile(meta.getAttribute(REGEX_KEY).get());
				}
			else
				{
				this.acceptRegex=null;
				}
			this.extend= meta.getIntAttribute(EXTEND_KEY).orElse(0);
			}
		else
			{
			this.fileHeader = null;
			this.acceptRegex = null;
			this.extend = 0;
			}
		}

	@Override
	public void fillHeader(final VCFHeader header) {
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		if(dict!=null && !SequenceDictionaryUtils.isGRCh38(dict)) {
			throw new IllegalArgumentException("VCF Sequence dictionary doesn't look like grch38 and database grch38 in "+getUri());
			}
		final String cmdHelp = " for regulomedb "+this.getUri() + " ranking-regex:"+(this.acceptRegex==null?"none":this.acceptRegex.pattern())+
		". The scoring scheme refers to the following supporting evidence for that particular location or variant id. "
		+ "In general, if more supporting data is available, the higher is its likelihood of being functional and hence receives a higher score (with 1 being higher and 7 being lower score)."
			;
		this.hdrRegulomeDbMeanScore = new VCFInfoHeaderLine(
				"REGDB_MEAN_SCORE",
				1,
				VCFHeaderLineType.Float,
				"Mean probability_score"+cmdHelp);
		this.hdrRegulomeDbMaxScore = new VCFInfoHeaderLine(
				"REGDB_MAX_SCORE",
				1,
				VCFHeaderLineType.Float,
				"Max probability_score "+cmdHelp);

		header.addMetaDataLine(this.hdrRegulomeDbMaxScore);
		header.addMetaDataLine(this.hdrRegulomeDbMeanScore);
		}
	
	@Override
	public List<VariantContext> annotate(VariantContext ctx) throws IOException {
		if(!isValid() || !hasContig(ctx)) return Collections.singletonList(ctx);
		final int start=Math.max(0,ctx.getStart()-this.extend);
		final int end=ctx.getEnd()+this.extend;

		TabixReader.Iterator r = super.tabixReader.query(contig(ctx),Math.max(start-1,0), end);
		double max_probability_score = 0.0;
		double sum_probability_score = 0.0;
		int count_probability_score = 0;
		for(;;) {
			String line = r.next();
			if(line==null) break;
			final Map<String, String> curr = this.fileHeader.toMap(CharSplitter.TAB.split(line));
			if(this.acceptRegex!=null && 
					   !this.acceptRegex.matcher(curr.get("ranking")).matches()
					   )
					{
					continue;
					}
			final int x0 = Integer.parseInt(curr.get("start")+1);
			final int x1 = Integer.parseInt(curr.get("end"));
			if(!CoordMath.overlaps(start, end, x0, x1)) continue;
			
			final double probability_score = Double.parseDouble(curr.get("probability_score"));
			max_probability_score = Math.max(probability_score, max_probability_score);
			sum_probability_score+=probability_score;
			count_probability_score++;
			}
		if(count_probability_score>0) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.attribute(this.hdrRegulomeDbMaxScore.getID(), max_probability_score);
			vcb.attribute(this.hdrRegulomeDbMeanScore.getID(), sum_probability_score/count_probability_score);
			return Collections.singletonList(vcb.make());
			}
		else
			{
			return Collections.singletonList(ctx);
			}
		
    	}

}
