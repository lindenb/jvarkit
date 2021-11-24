/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.roh;


import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
/*
BEGIN_DOC

 
## Example


```bash
java -jar dist/vcfpar.jar   in.vcf > out.vcf
```

END_DOC
 */
@Program(name="vcfroh",
	description="VCF ROH",
	keywords={"vcf","roh"},
	creationDate="20211123",
	modificationDate="20211123"
	)
public class VcfROH extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfROH.class).make();

	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;
	@Parameter(names={"--score"},description="HOM_REF/HOM_VAR score")
	private double hom_score = 0.025;
	@Parameter(names={"--nc2hr"},description="treat NO_CALL to HOM_REF")
	private boolean nocall_to_homref = false;

	private class Sample {
		final String name;
		final int index;
		String prev_contig = null;
		int start = -1;
		int end  = -1;
		double score=0;
		int count = 0;
		Sample(final String sn,int index) {
			this.name= sn;
			this.index = index;
			}
		private void dump(final PrintWriter w) {
			if(StringUtils.isBlank(this.prev_contig)) return;
			if(this.start>=this.end) return;
			w.print(this.prev_contig);
			w.print("\t");
			w.print(this.start);
			w.print("\t");
			w.print(this.end);
			w.print("\t");
			w.print(CoordMath.getLength(this.start, this.end));
			w.print("\t");
			w.print(this.name);
			w.print("\t");
			w.print(this.count);
			w.print("\t");
			w.print(this.score);
			w.println();
		}
		
		void visit(final PrintWriter w,final VariantContext ctx) {
			final Genotype gt = ctx.getGenotype(this.index);
			final double dx;
			switch(gt.getType()) {
				case HOM_VAR: //continue
				case HOM_REF: dx = VcfROH.this.hom_score; break;
				case NO_CALL: {
					if (VcfROH.this.nocall_to_homref) {
						dx = VcfROH.this.hom_score;
						break;
						}
					//continue
					}
				default: dx = 1.0 - VcfROH.this.hom_score; break;
				}
			if(!ctx.getContig().equals(this.prev_contig) || this.score+dx <= 0) {
				dump(w);
				this.prev_contig = ctx.getContig();
				this.start = ctx.getStart();
				this.score=0;
				this.count  = 0;
				}
			this.end = ctx.getEnd();
			this.score += dx;
			this.count++;
			}
		void finish(PrintWriter w) {
			dump(w);
		}
	}
	
	public VcfROH() {
		}
	
	
	@Override
	public int doWork(List<String> args) {
		try {
			final String input = oneAndOnlyOneFile(args);
			try(VCFIterator iter = input==null?new VCFIteratorBuilder().open(stdin()):new VCFIteratorBuilder().open(input)) {
				final VCFHeader header= iter.getHeader();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				final OrderChecker<VariantContext> checker = new OrderChecker<VariantContext>(dict,false);
				final List<Sample> samples = new ArrayList<>(header.getNGenotypeSamples());
				for(String sn:header.getSampleNamesInOrder()) {
					final Sample sample =new Sample(sn, header.getSampleNameToOffset().get(sn));
					samples.add(sample);
				}
				
				try(PrintWriter w =  super.openPathOrStdoutAsPrintWriter(this.output)) {
					while(iter.hasNext()) {
						final VariantContext ctx = checker.apply(iter.next());
						for(Sample sample:samples) sample.visit(w,ctx);
						}
					for(Sample sample:samples) sample.finish(w);
					w.flush();
					}

				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		
		}
	


	public static void main(final String[] args)
		{
		new VcfROH().instanceMainWithExit(args);
		}
	
	}
