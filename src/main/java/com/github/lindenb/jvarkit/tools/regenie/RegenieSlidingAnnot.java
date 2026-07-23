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
package com.github.lindenb.jvarkit.tools.regenie;

import java.io.PrintWriter;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

 The aim of this  class is to produce a file for regenie containing sliding windows with the following header:
 
<pre>"CONTIG","POS","ID","GENE","ANNOTATION","SCORE","CADD","FREQ","SINGLETON"</pre>
 

## Example


END_DOC
**/
@Program(name="regenieslidingannot",
description="Create annotation files for regenie using sliding annotations",
keywords={"vcf","regenie","burden"},
creationDate="20250311",
modificationDate="20250320",
jvarkit_amalgamion = true,
generate_doc = true
)
public class RegenieSlidingAnnot extends AbstractRegenieAnnot {
	private static final Logger LOG = Logger.of(RegenieSlidingAnnot.class);
	@Parameter(names = {"--window-size","-w"}, description = "window size. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter=DistanceParser.StringConverter.class,required =true)
	private int window_size=-1;
	@Parameter(names = {"--window-shift","-s"}, description = "window shift. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter=DistanceParser.StringConverter.class,required = true)
	private int window_shift=-1;

	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected VCFHeader initVcfHeader(final VCFHeader h) {
		if(window_shift<1 || window_shift>window_size) throw new IllegalArgumentException("shift>size");
		return super.initVcfHeader(h);
		}

	private String getPredictionName() {
		return "sliding_"+window_size;
	}

	@Override
	protected void dump(final PrintWriter w,final VariantContext ctx) throws Exception {
		final byte is_singleton = isSingleton(ctx);
		final double freq = getFrequency(ctx);
		int  win_pos = 1 + (((int)(ctx.getStart()/(double)this.window_size)) * this.window_size);
		do {
			final Variation v = new Variation();
			v.contig = fixContig(ctx.getContig());
			v.pos = ctx.getStart();
			v.id = makeID(ctx);
			v.gene = ctx.getContig()+ "_" + (win_pos) + "_" + (win_pos - 1 + this.window_size) ;
			v.prediction = getPredictionName();
			v.score =(isIgnoringMaskScore()?OptionalDouble.of(1.0): OptionalDouble.empty());
			v.cadd = getCaddScore(ctx);
			v.is_singleton = is_singleton;
			v.frequency = freq;
			print(w,v);
			win_pos+=this.window_shift;
			} while(CoordMath.overlaps(win_pos, win_pos+this.window_size-1, ctx.getStart(), ctx.getEnd()));
		}

	public static void main(final String[] args) {
		new RegenieSlidingAnnot().instanceMainWithExit(args);
	}
}
