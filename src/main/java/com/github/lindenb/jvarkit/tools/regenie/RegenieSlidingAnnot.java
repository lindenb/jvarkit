package com.github.lindenb.jvarkit.tools.regenie;

import java.io.PrintWriter;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC


END_DOC
**/
@Program(name="regenieslidingannot",
description="Create annotation files for regenie using sliding annotations",
keywords={"vcf","regenie","burden"},
creationDate="20250311",
modificationDate="20250311",
generate_doc = false
)
public class RegenieSlidingAnnot extends AbstractRegenieAnnot {
	private static final Logger LOG = Logger.build(RegenieSlidingAnnot.class).make();
	@Parameter(names = {"--window-size"}, description = "window size. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter=DistanceParser.StringConverter.class,required =true)
	private int window_size=-1;
	@Parameter(names = {"--window-shift"}, description = "window shift. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter=DistanceParser.StringConverter.class,required = true)
	private int window_shift=-1;

	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected VCFHeader initVcfHeader(VCFHeader h) {
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
			v.score = OptionalDouble.empty();
			v.cadd = getCaddScore(ctx);
			v.is_singleton = is_singleton;
			v.frequency = freq;
			print(w,v);
			win_pos+=this.window_shift;
			} while(CoordMath.overlaps(win_pos, win_pos+this.window_size-1, ctx.getStart(), ctx.getEnd()));
		}

	public static void main(String[] args) {
		new RegenieSlidingAnnot().instanceMainWithExit(args);
	}
}
