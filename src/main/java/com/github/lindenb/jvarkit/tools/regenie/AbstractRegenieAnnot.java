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
package com.github.lindenb.jvarkit.tools.regenie;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
 * Base class for Regenie** Generator
 *
 */
public abstract class AbstractRegenieAnnot extends Launcher {
	private static final String CADD_PHRED = "CADD_PHRED";
	private static final String GNOMAD_AF = "gnomad_genome_AF_NFE";
	@Parameter(names = "-o", description = OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;

	@Parameter(names = "-f", description = "comma separated of Allele frequencies , I will use the highest to discard frequent variants.")
	private String freqStr="0.01";

	protected static class Variation {
		String contig;
		int pos;
		String id;
		String gene;
		String prediction;
		OptionalDouble score=OptionalDouble.empty();
		OptionalDouble cadd=OptionalDouble.empty();
		int is_singleton=0;
		double frequency=0.0;
	}
	

	protected byte isSingleton(final VariantContext ctx) {
		if(!ctx.hasGenotypes()) throw new IllegalArgumentException("GT missing");
		Genotype single=null;
		for(Genotype g:ctx.getGenotypes()) {
			if(!g.hasAltAllele()) continue;
			if(single!=null) return 0;
			single =g;
		}
		return (byte)(single==null?0:1);
		}
	
	protected double getFrequency(final VariantContext ctx) {
		if(!ctx.hasGenotypes()) throw new IllegalArgumentException("GT missing");
		int AC=0;
		double AN=0;
		final Allele the_alt = ctx.getAlternateAllele(0);
		for(Genotype gt: ctx.getGenotypes()) {
			for(Allele alt: gt.getAlleles()) {
				if(the_alt.equals(alt)) {
					AC++;
					}
				AN++;//consider NO_CALL=REF
				}
			}
		double freq = AC/AN;
		if(ctx.hasAttribute(GNOMAD_AF)) {
			freq = Math.max(ctx.getAttributeAsDouble(GNOMAD_AF, 0.0), freq);
			}
		return freq;
		}
	
	protected abstract void dump(final PrintWriter w,final VariantContext ctx) throws Exception;
		
	
	protected void print(final PrintWriter w,final Variation ctx)  {
		w.print(ctx.contig);
		w.print("\t");
		w.print(ctx.pos);
		w.print("\t");
		w.print(ctx.id);
		w.print("\t");
		w.print(ctx.gene);
		w.print("\t");
		w.print(ctx.prediction);
		w.print("\t");
		w.print(ctx.score.orElse(1.0));
		w.print("\t");
		w.print(ctx.cadd.orElse(0.0));
		w.print("\t");
		w.print(ctx.frequency);
		w.print("\t");
		w.print(ctx.is_singleton);
		w.println();
	}
	

	protected String makeID(final VariantContext vc) {
		return String.join(":",
			fixContig(vc.getContig()),
			String.valueOf(vc.getStart()),
			vc.getReference().getDisplayString(),
			vc.getAlternateAllele(0).getDisplayString()
			);
		}
	
	protected String fixContig(final String ctg) {
		//if(ctg.equals("X") || ctg.equals("chrX")) return "23";
		//if(ctg.equals("Y") || ctg.equals("chrY")) return "24";
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	private boolean keepVariant(final double max_freq,final VariantContext ctx) {
		return getFrequency(ctx) <= max_freq;
		}

	protected OptionalDouble getCaddScore(final VariantContext ctx) {
		if (ctx.hasAttribute(CADD_PHRED)) {
				final String s = ctx.getAttributeAsString(CADD_PHRED, ".");
				if (!(s.equals(".") || StringUtils.isBlank(s))) {
					return OptionalDouble.of(Double.valueOf(s));
				}
			}
		return OptionalDouble.empty();
		}
	
	protected VCFHeader initVcfHeader(final VCFHeader h) {
		return h;
	}
	
	protected abstract Logger getLogger();
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final double freq = Arrays.stream(CharSplitter.COMMA.split(this.freqStr.trim())).mapToDouble(S->Double.parseDouble(S)).max().orElse(0.01);
			final String input = oneFileOrNull(args);
			
			try (VCFIterator iter = VCFUtils.createVCFIterator(input)) {
				initVcfHeader(iter.getHeader());
				try (PrintWriter w = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					final String[] header_line=new String[]{"CONTIG","POS","ID","GENE","ANNOTATION","SCORE","CADD","FREQ","SINGLETON"};
					w.println(String.join("\t", header_line));
					
					
					while (iter.hasNext()) {
						final VariantContext vc = iter.next();
						if (vc.getNAlleles() != 2)
							throw new IOException(vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles());
						if(!keepVariant(freq,vc)) continue;
						dump(w, vc);
					} // end while
				w.flush();
				}
			}

			return 0;
		} catch (final Throwable err) {
			getLogger().error(err);
			return -1;
		}
	}

}
