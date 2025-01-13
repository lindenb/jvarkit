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
package com.github.lindenb.jvarkit.tools.prs;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**

BEGIN_DOC

# buid tabix database

```
$ head AR_test.txt | column -t
#CHROM  POS        ID          REF  ALT  EFFECT_ALL  EFFECT  
chr1    5550460   rs249409   G    A    G           0.052
chr1    10918306  rs69301    G    T    T           0.15
chr2    2126390   rs13617   G    A    A           0.1
chr2    4407256   rs476   G    T    G           0.071
chr6    1605860  rs1564348   T    C    T           0.014
chr6    2093141   rs18   G    A    G           0.057

LC_ALL=C sort -T . -t $'\t' -k1,1 -k2,2n  AR_test.txt | bgzip > AR_test.txt.gz
tabix -s1 -b 2 -e 2  -c '#' 20220915_AR_test.txt.gz
```

and then

```
bcftools view in.bcf | java -jar dist/vcfsamplesprs.jar -S AR_test.txt.gz
```

END_DOC
*/
@Program(
		name="vcfsamplesprs",
		description="another program for @AntoineRimbert",
		keywords={"vcf","indel"},
		creationDate = "20220915",
		modificationDate = "20220915",
		generate_doc = true
		)
public class VcfSamplesPRS extends Launcher {
	private static final Logger LOG = Logger.build(VcfSamplesPRS.class).make();
	private static final boolean IGNORE_REF_STATE = true;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	
	@Parameter(names={"--score","--scores","-S"},description="tabix indexed scores. CHROM(tab)POS(tab)ID(tab)REF(tab)ALT(tab)EFFECT_ALL(tab)EFFECT  ",required = true)
	private String scoreVCFfile = null;

	
	private static class Sample {
		final int index;
		final String name;
		double score=0;
		long n_nocall = 0L;
		Sample(int index,String name) {
			this.index = index;
			this.name = name;
			}
		}
	
	private static class Score {
		final String contig;
		final int pos;
		final String id;
		final Allele REF;
		final Allele ALT;
		final Allele EFFECT_ALLELE;
		final double score;
		Score(final String[] tokens) {
			if(tokens.length!=7) throw new IllegalArgumentException("expected 7 tokens in "+String.join("(tab)", tokens));
			contig = tokens[0];
			pos = Integer.parseInt(tokens[1]);
			id =  tokens[2];
			REF =  Allele.create(tokens[3],true);
			ALT =  Allele.create(tokens[4],false);
			if(ALT.equals(REF,IGNORE_REF_STATE)) {
				throw new IllegalArgumentException("REF==ALT in "+String.join("(tab)", tokens));
				}

			EFFECT_ALLELE =  Allele.create(tokens[5],false);
			if(!(EFFECT_ALLELE.equals(REF,IGNORE_REF_STATE) || EFFECT_ALLELE.equals(ALT,IGNORE_REF_STATE))) {
				throw new IllegalArgumentException("EFFECT_ALL is not in REF or ALT in "+String.join("(tab)", tokens));
				}
			score = Double.parseDouble(tokens[6]);
			}
		@Override
		public String toString() {
			return contig+":"+pos+":"+id+":"+REF.getDisplayString()+"/"+ALT.getDisplayString();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			try(TabixReader tabix = new TabixReader(this.scoreVCFfile)) {
				final ContigNameConverter ctgConverter = ContigNameConverter.fromContigSet(tabix.getChromosomes());
				try (VCFIterator r = super.openVCFIterator(input)){
					final VCFHeader h = r.getHeader();
					final List<Sample> samples = new ArrayList<>(h.getNGenotypeSamples());
					h.getSampleNameToOffset().entrySet().stream().forEach(KV->{
						samples.add(new Sample(
							KV.getValue().intValue(),
							KV.getKey()
							));
						});
					while(r.hasNext()) {
						final VariantContext ctx = r.next();
						if(ctx.getNAlleles()!=2) {
							LOG.warn("skipping multiallelic "+ctx.getContig()+":"+ctx.getStart());
							continue;
							}
						final String ctgTabix = ctgConverter.apply(ctx.getContig());
						if(StringUtils.isBlank(ctgTabix)) {
							LOG.warn("chromosome not in "+this.scoreVCFfile+" "+ctx.getContig()+":"+ctx.getStart());
							continue;
							}
						final List<Score> scores = new ArrayList<>();
						TabixReader.Iterator iter =  tabix.query(ctgTabix,Math.max(0,ctx.getStart()-1), ctx.getEnd()+1);
						for(;;) {
							
							final String line = iter.next();
							if(line==null) break;
							final Score score = new Score(CharSplitter.TAB.split(line));
							if(score.pos!=ctx.getStart()) continue;
							if(score.REF.equals(ctx.getAlleles().get(1),IGNORE_REF_STATE) &&
								score.ALT.equals(ctx.getReference(),IGNORE_REF_STATE)) {
									LOG.warn("score ALT/REF swapped for "+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString());
									continue;
									}
							
							if(!score.REF.equals(ctx.getReference(),IGNORE_REF_STATE)) continue;
							if(!score.ALT.equals(ctx.getAlleles().get(1),IGNORE_REF_STATE)) continue;
							scores.add(score);
							}
						if(scores.isEmpty()) {
							LOG.warn("no score for "+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString());
							continue;
							}
						if(scores.size()!=1) {
							LOG.warn("skipping multiple scores "+scores);
							continue;
							}
						LOG.warn("got "+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString());
						final Score score = scores.get(0);
						
						for(Sample sn:samples) {
							final Genotype g = ctx.getGenotype(sn.index);
							switch(g.getType()) {
								case HET:
									sn.score += score.score;
									break;
								case HOM_VAR:
									if(!score.REF.equals(score.EFFECT_ALLELE,IGNORE_REF_STATE))
										{
										sn.score+= 2.0 * score.score; 
										}
									break;
								case HOM_REF:{
									if(score.REF.equals(score.EFFECT_ALLELE,IGNORE_REF_STATE)) {
										sn.score+= 2.0 * score.score; 
										}
									break;
									}
								case MIXED:
								case NO_CALL:
								case UNAVAILABLE:
									sn.n_nocall++;
									break;
								}
							}
						}
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					pw.println("#SAMPLE\tSCORE\tNO_CALLS");
					Collections.sort(samples,(A,B)->A.name.compareTo(B.name));
					for(Sample sn:samples) {
						pw.print(sn.name);
						pw.print("\t");
						pw.print(sn.score);
						pw.print("\t");
						pw.print(sn.n_nocall);
						pw.println();
						}
					pw.flush();
					}
				}
			}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	
	public static void main(String[] args) {
		new VcfSamplesPRS().instanceMainWithExit(args);
	}

}
