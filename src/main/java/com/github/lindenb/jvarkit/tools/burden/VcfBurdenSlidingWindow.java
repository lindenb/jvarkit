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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.iterator.SlidingWindowIterator;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

# Example

```
$ java -jar dist/vcfburdenslidingwindow.jar --pedigree ./src/test/resources/test_vcf01.ped -t 1 ./src/test/resources/test_vcf01.vcf  | head

#chrom	start0	end	name	length	p-value	affected_alt	affected_hom	unaffected_alt	unaffected_hom	variants.count
1	832199	833200	1:832200-833200	1001	1.0	0	3	1	2	1
1	832499	833500	1:832500-833500	1001	1.0	0	3	1	2	1
1	832799	833800	1:832800-833800	1001	1.0	0	3	1	2	1
1	839999	841000	1:840000-841000	1001	1.0	0	3	1	2	1
1	840299	841300	1:840300-841300	1001	1.0	0	3	1	2	1
1	840599	841600	1:840600-841600	1001	1.0	0	3	1	2	1
1	849299	850300	1:849300-850300	1001	1.0	0	3	1	2	1
1	849599	850600	1:849600-850600	1001	1.0	1	2	1	2	2
1	849899	850900	1:849900-850900	1001	1.0	1	2	1	2	2
```

END_DOC
*/

@Program(name="vcfburdenslidingwindow",
description="Run Burden Sliding window",
keywords={"vcf","burden","case","control"},
creationDate="20190920",
modificationDate="20190920"
)
public class VcfBurdenSlidingWindow
extends AbstractVcfBurden
{
	private static final Logger LOG = Logger.build(VcfBurdenSlidingWindow.class).make();
	@Parameter(names={"-f","--filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-t","--treshold"},description="fisher-test treshold. Discard results greater than this value.")
	private double fisherTreshold = 1e-5;
	@Parameter(names={"-w","--window-size"},description="Window size",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int window_size=1_000;
	@Parameter(names={"-s","--window-shift"},description="Window shift",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int window_shift=300;
	@Parameter(names={"-C","--contig"},description="limit to this contig")
	private String limitContig = null;


	private boolean accept(final VariantContext ctx) {
		if(!variantFilter.test(ctx))return false;
		return true;
	}

	@Override
	protected void runBurden(final PrintWriter pw, final VCFReader vcfReader, final VariantContextWriter vcw) throws IOException {
		if(this.window_size<=0) throw new IllegalArgumentException("bad window size: "+this.window_size);
		if(this.window_shift<=0) throw new IllegalArgumentException("bad window shift:" + this.window_shift);
		
		final SAMSequenceDictionary vcfDict = SequenceDictionaryUtils.extractRequired(vcfReader.getHeader());
		
		pw.print("#chrom");
		pw.print("\t");
		pw.print("start0");
		pw.print("\t");
		pw.print("end");
		pw.print("\t");
		pw.print("name");
		pw.print("\t");
		pw.print("length");
		pw.print("\t");
		pw.print("p-value");
		pw.print("\t");
		pw.print("affected_alt");
		pw.print("\t");
		pw.print("affected_hom");
		pw.print("\t");
		pw.print("unaffected_alt");
		pw.print("\t");
		pw.print("unaffected_hom");
		pw.print("\t");
		pw.print("variants.count");
		pw.println();

		
		final ProgressFactory.Watcher<Locatable> progress = ProgressFactory.newInstance().logger(LOG).dictionary(vcfDict).build();
		
		
		
		CloseableIterator<VariantContext> iter;
		if(!StringUtils.isBlank(this.limitContig)) {
			final SAMSequenceRecord ssr= vcfDict.getSequence(this.limitContig);
			if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(limitContig, vcfDict);
			iter = vcfReader.query(ssr);
			}
		else
			{	
			iter = vcfReader.iterator();
			}
	
		
		final SlidingWindowIterator<VariantContext> iter2 = new SlidingWindowIterator<>(
				iter.stream().filter(CTX->accept(CTX)).iterator(),
				this.window_size,
				this.window_shift
				);
		while(iter2.hasNext()) {
			final Map.Entry<? extends Locatable, List<VariantContext>> entry = iter2.next();
			
			final Locatable W =  progress.apply(entry.getKey());
			final FisherResult fisher = runFisher(entry.getValue());
			if(fisher.p_value> this.fisherTreshold) continue;
			
			pw.print(W.getContig());
			pw.print("\t");
			pw.print(W.getStart()-1);
			pw.print("\t");
			pw.print(W.getEnd());
			pw.print("\t");
			pw.print(new SimpleInterval(W).toString());
			pw.print("\t");
			pw.print(W.getLengthOnReference());
			pw.print("\t");
			pw.print(fisher.p_value);
			pw.print("\t");
			pw.print(fisher.affected_alt);
			pw.print("\t");
			pw.print(fisher.affected_hom);
			pw.print("\t");
			pw.print(fisher.unaffected_alt);
			pw.print("\t");
			pw.print(fisher.unaffected_hom);
			pw.print("\t");
			pw.print(entry.getValue().size());
			pw.println();

			
			
			if(vcw!=null) {
				for(final VariantContext ctx:entry.getValue()) {
					vcw.add(new VariantContextBuilder(ctx).
							attribute(BURDEN_KEY, VCFUtils.escapeInfoField(W.toString())).
							make());
				}
			}
		
		}
		
	
	iter.close();
	progress.close();	
	}
	
public static void main(final String[] args) {
	new VcfBurdenSlidingWindow().instanceMainWithExit(args);
	}
}
