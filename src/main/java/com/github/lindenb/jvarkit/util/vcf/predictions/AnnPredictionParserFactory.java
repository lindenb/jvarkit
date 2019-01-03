/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.vcf.predictions;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class AnnPredictionParserFactory extends AbstractPredictionParserFactory<AnnPredictionParser,AnnPredictionParserFactory>
	{
	public AnnPredictionParserFactory() {
		super(AnnPredictionParser.getDefaultTag());
		}
	public AnnPredictionParserFactory(final VCFHeader header) {
		super(AnnPredictionParser.getDefaultTag(),header);
		}
	
	@Override
	public AnnPredictionParser get() {
		return new AnnPredictionParser(getHeader(), getTag());
	}
	
	/** create a parser without the VCF header */
	public AnnPredictionParser createDefaultParser() 
		{
		final VCFHeader tmpheader=new VCFHeader();
		final VCFInfoHeaderLine info=new VCFInfoHeaderLine(
				AnnPredictionParser.getDefaultTag(),
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				""
				);
		tmpheader.addMetaDataLine(info);
		return new AnnPredictionParser(tmpheader, getTag());
		}
}
