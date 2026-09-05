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
package com.github.lindenb.jvarkit.tools.vcfcoocurrence;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.LogicalOp;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherCasesControls;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.CasesControls;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
/**
BEGIN_DOC

Burden test for pair of variants

```
$ java -jar dist/jvarkit.jar vcfcoocurrence --cases S1,S2 --controls S4,S5 src/test/resources/rotavirus_rf.vcf.gz
#CHROM1 START1  REF1    ALT2    CHROM2  START1  REF2    ALT2    CASE_ALT        CASE_REF        CTRL_ALT        CTRL_REF        PVALUE
RF01    970     A       C       RF03    1688    T       G       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF03    1708    G       T       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF04    1900    A       C       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF06    668     A       G       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF07    952     A       G       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF08    926     A       C       0       2       1       1       0.9999999999999999
RF01    970     A       C       RF09    317     C       A       0       2       1       1       0.9999999999999999
RF02    251     A       T       RF02    1726    T       G       1       1       0       2       0.9999999999999999
```

END_DOC
 */
// https://github.com/hail-is/hail/blob/1a861505c1fc2ea3c9d7b32a47be7af10d13907c/hail/src/main/scala/is/hail/experimental/package.scala#L38-L84
// https://github.com/broadinstitute/gnomad-browser/blob/main/graphql-api/src/queries/variant-cooccurrence-queries.ts#L85
@Program(name="vcfcoocurrence",
	description="fisher-test of coocurrence of variants",
	keywords={"vcf","fisher","burden"},
	jvarkit_amalgamion =  true,
	creationDate = "2026-09-04",
	modificationDate = "2026-09-05"
	)
public class VcfCoOccurrence extends Launcher {
	private static final Logger LOG = Logger.of(VcfCoOccurrence.class);
	@Parameter(names={"--treshold","-t"},description="fisher treshold, don't print pairs with fisher test > t")
	private double treshold  = 1.0;
	@Parameter(names={"--internal","-i"},description="internal fisher treshold, don't use variant where fisher exact test is lower than 'x'")
	private double internal_treshold  = 0.0;
	@Parameter(names={"--operator","--logical","--op"},description="logical operator to combine to variants. e.g: AND= both genotypes must be carried by the same sample. OR: one of the genotype must be carried by the same sample, etc...")
	private LogicalOp logical_op = LogicalOp.AND;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private CasesControls caseControls = new CasesControls();

	private static class Variant implements Locatable {
		final String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		final BitSet has_alt;
		double internal_fisher=0.0;
		Variant(final VariantContext ctx,final Allele alt) {
			this.contig = ctx.getContig();
			this.pos = ctx.getStart();
			this.ref = ctx.getReference();
			this.alt = alt;
			this.has_alt = new BitSet(ctx.getNSamples());
			for(int i=0;i< ctx.getNSamples();i++) {
				if(ctx.getGenotype(i).getAlleles().contains(alt)) {
					this.has_alt.set(i);
					}
				}
			this.internal_fisher = 0;
			}
		
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return pos;
			}
		@Override
		public int getEnd() {
			return pos + ref.length()-1;
			}
		}
	

	private List<Variant> toVariants(final VariantContext ctx,
			final int[] cases_indexes,
			final int[] ctrls_indexes
			) {
		final int[] counts = new int[4];
		final List<Variant> list =new ArrayList<Variant>(ctx.getNAlleles());
		for(Allele alt: ctx.getAlternateAlleles()) {
			if(alt.equals(Allele.SPAN_DEL)) continue;
			if(alt.equals(Allele.NO_CALL)) continue;
			final Variant v = new Variant(ctx,alt);
			if(v.has_alt.nextSetBit(0)==-1) continue;//all bits are off
			
			Arrays.fill(counts, 0);
			for(int x: cases_indexes) {
				if(v.has_alt.get(x)) {
					counts[0]++;
					}
				else
					{
					counts[1]++;
					}
				}
			for(int x: ctrls_indexes) {
				if(v.has_alt.get(x)) {
					counts[2]++;
					}
				else
					{
					counts[3]++;
					}
				}
			if(counts[0]==0 && counts[2]==0) continue; //no alt
			final FisherExactTest fisher = FisherExactTest.compute(counts);
			v.internal_fisher = fisher.getAsDouble();
			if(v.internal_fisher < this.internal_treshold) continue;
			list.add(v);
			}
		return list;
		}
	

	private void scan(
			PrintWriter w,
			final VCFIterator iter,
			final int[] cases_indexes,
			final int[] ctrls_indexes
			) {
		final int[] counts = new int[4];
		final List<Variant> variants = iter.stream()
				.flatMap(V->toVariants(V,cases_indexes,ctrls_indexes).stream())
				.collect(Collectors.toList());
		
		for(int i=0;i+1 < variants.size();i++) {
			final Variant vi =  variants.get(i);
			for(int j=i+1;j < variants.size();j++) {
				final Variant vj =  variants.get(j);
				Arrays.fill(counts, 0);
				for(int x: cases_indexes) {
					if(this.logical_op.test( vi.has_alt.get(x) , vj.has_alt.get(x))) {
						counts[0]++;
						}
					else
						{
						counts[1]++;
						}
					}
				for(int x: ctrls_indexes) {
					if(this.logical_op.test(vi.has_alt.get(x) , vj.has_alt.get(x))) {
						counts[2]++;
						}
					else
						{
						counts[3]++;
						}
					}
				if(counts[0]==0 && counts[2]==0) continue;
				final FisherExactTest fisher = FisherExactTest.compute(counts);
				final double pvalue= fisher.getAsDouble();
				if(pvalue > this.treshold) continue;
				w.print(vi.getContig());
				w.print("\t");
				w.print(vi.getStart());
				w.print("\t");
				w.print(vi.ref.getDisplayString());
				w.print("\t");
				w.print(vi.alt.getDisplayString());
				w.print("\t");
				w.print(vi.internal_fisher);
				w.print("\t");
				w.print(vj.getContig());
				w.print("\t");
				w.print(vj.getStart());
				w.print("\t");
				w.print(vj.ref.getDisplayString());
				w.print("\t");
				w.print(vj.alt.getDisplayString());
				w.print("\t");
				w.print(vj.internal_fisher);
				w.print("\t");
				w.print(counts[0]);
				w.print("\t");
				w.print(counts[1]);
				w.print("\t");
				w.print(counts[2]);
				w.print("\t");
				w.print(counts[3]);
				w.print("\t");
				w.print(pvalue);
				w.println();
				}
			}
		}


	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input= oneFileOrNull(args);
			
			try(VCFIterator iter= StringUtils.isBlank(input)?
					new VCFIteratorBuilder().open(stdin()):
					new VCFIteratorBuilder().open(input)
					)
				{
				final VCFHeader header = iter.getHeader();
				this.caseControls.load();
				this.caseControls.retain(header);
				this.caseControls.checkHaveCasesControls();
				final Map<String,Integer> sample2offset = header.getSampleNameToOffset();
				final int[] cases = this.caseControls.getCases().stream().map(S->sample2offset.get(S)).filter(I->I!=null).mapToInt(I->I.intValue()).sorted().toArray();
				final int[] controls =this. caseControls.getControls().stream().map(S->sample2offset.get(S)).filter(I->I!=null).mapToInt(I->I.intValue()).sorted().toArray();
				
				if(cases.length==0) {
					LOG.error("No cases in VCF");
					return -1;
					}
				if(controls.length==0) {
					LOG.error("No controls in VCF");
					return -1;
					}
				
				try(PrintWriter w =  super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					w.print("#CHROM1");
					w.print("\t");
					w.print("START1");
					w.print("\t");
					w.print("REF1");
					w.print("\t");
					w.print("ALT1");
					w.print("\t");
					w.print("FISHER1");
					w.print("\t");
					w.print("CHROM2");
					w.print("\t");
					w.print("START2");
					w.print("\t");
					w.print("REF2");
					w.print("\t");
					w.print("ALT2");
					w.print("\t");
					w.print("FISHER2");
					w.print("\t");
					w.print("CASE_ALT");
					w.print("\t");
					w.print("CASE_REF");
					w.print("\t");
					w.print("CTRL_ALT");
					w.print("\t");
					w.print("CTRL_REF");
					w.print("\t");
					w.print("FISHER");
					w.println();
					
					scan(
						w,
						iter,
						cases,
						controls
						);
					w.flush();
					}
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(String[] args) {
		new VcfCoOccurrence().instanceMainWithExit(args);
	}
}
