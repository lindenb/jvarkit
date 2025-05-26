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
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.math.MinMaxDouble;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 *  Variant Annotator using BigWig  files
 */
public class BigWigVariantAnnotator implements VariantAnnotator {
	private final BigWigReader bigWigReader;
	private String tag;
	private VCFInfoHeaderLine averageInfo = null;
	private VCFInfoHeaderLine minInfo = null;
	private VCFInfoHeaderLine maxInfo = null;
	private long count_variants =0L;//for System.gc call

	
	public BigWigVariantAnnotator(final String bigwigUri)  throws IOException {
		this.bigWigReader = new BigWigReader(bigwigUri);
		
		this.tag = bigwigUri;
		int i=this.tag.lastIndexOf(File.separator);
		if(i!=-1) this.tag=this.tag.substring(i+1);
		i=this.tag.lastIndexOf('.');
		if(i!=-1) this.tag=this.tag.substring(0,i);
		if(StringUtils.isBlank(tag)) tag = "TAG";
		}
	
	public String getTag() {
		return tag;
		}
	
	public void setTag(final String tag) {
		this.tag = tag;
		if(StringUtils.isBlank(tag)) throw new IllegalArgumentException("tag shouldn't be empty");
		}
	
	@Override
	public void fillHeader(final VCFHeader header) {
		final String tag = getTag();
		if(StringUtils.isBlank(tag)) throw new IllegalArgumentException("tag shouldn't be empty");
		
		
		this.averageInfo = new VCFInfoHeaderLine(
				tag,1,
				VCFHeaderLineType.Float,
				"Values from bigwig file: "+ this.bigWigReader.getURI()+ ". Average value if the variant covers multiple positions"
				);
		header.addMetaDataLine(this.averageInfo);
		
		this.minInfo = new VCFInfoHeaderLine(
				tag+"_MIN",1,
				VCFHeaderLineType.Float,
				"Min Value from bigwig file: "+ this.bigWigReader.getURI()+ " if the variant covers multiple positions and min!=average"
				);
		header.addMetaDataLine(this.minInfo);

		this.maxInfo = new VCFInfoHeaderLine(
				tag+"_MAX",1,
				VCFHeaderLineType.Float,
				"Min Value from bigwig file: "+ this.bigWigReader.getURI()+ " if the variant covers multiple positions and max!=average"
				);
		header.addMetaDataLine(this.maxInfo);
		}

	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
		final Locatable loc;
		// https://github.com/lindenb/jvarkit/issues/165
		if(ctx.isIndel()) {
			loc = new SimpleInterval(ctx.getContig(),ctx.getStart()+1,Math.max(ctx.getStart()+1,ctx.getEnd()));
			}
		else
			{
			loc = ctx;
			}

		
		final MinMaxDouble minMax= new MinMaxDouble();
		double sum=0.0;
		try(final CloseableIterator<BigWigReader.WigItem> iter=this.bigWigReader.query(loc)) {
			while(iter.hasNext())
				{
				final BigWigReader.WigItem item=iter.next();
				if(!item.overlaps(loc)) continue;//paranoid
				final float v=item.getValue();
				minMax.accept(v);
				sum+=v;
				}
			}
		if(minMax.isEmpty()) return Collections.singletonList(ctx);
		final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
		final double average = sum/minMax.getCount();
		vcb.attribute(this.averageInfo.getID(), average);
		if(minMax.getCount()>1L) {
			if(average!=minMax.getMinAsDouble()) vcb.attribute(this.minInfo.getID(), minMax.getMinAsDouble());
			if(average!=minMax.getMaxAsDouble()) vcb.attribute(this.maxInfo.getID(), minMax.getMaxAsDouble());
			}
		if(++count_variants%1000L==0)
			{
			System.gc();
			}
		
		return Collections.singletonList(vcb.make());
	}

	@Override
	public void close() {
		this.bigWigReader.close();
	}

}
