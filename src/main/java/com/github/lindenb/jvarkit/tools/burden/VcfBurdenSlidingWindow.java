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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.FilteringVariantContextIterator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

# Motivation

apply fisher test on VCF using a sliding window

# Example

```
$ java -jar dist/jvarkit.jar vcfburdenslidingwindow --cases cases.txt --controls controls.txt -t 1 ./src/test/resources/test_vcf01.vcf  | head

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
description="apply fisher test on VCF using a sliding window",
keywords={"vcf","burden","case","control"},
creationDate="20190920",
modificationDate="20231213",
jvarkit_amalgamion = true
)
public class VcfBurdenSlidingWindow extends Launcher {
	private static final Logger LOG = Logger.of(VcfBurdenSlidingWindow.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private CasesControls casesControls = new CasesControls();
	@Parameter(names={"-f","--filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,splitter = NoSplitter.class, converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-t","--treshold"},description="fisher-test treshold. Discard results greater than this value.",splitter=NoSplitter.class,converter=FractionConverter.class)
	private double fisherTreshold = 1.0;
	@Parameter(names={"-w","--window-size"},description="Window size",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int window_size=1_000;
	@Parameter(names={"-s","--window-shift"},description="Window shift",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int window_shift=300;
	@Parameter(names={"-C","--contig"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"-save-vcf","--save-vcf"},description="Save Matching variants for the best p-value into that VCF.")
	private Path outputVcfPath = null;

	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(window_size < 1 || window_shift < 1 || window_size< window_shift) {
				LOG.error("bad window size / window_shift");
				return -1;
				}
			final String input  = oneFileOrNull(args);

			double best_pvalue= 1.0;
			this.casesControls.load();
			this.casesControls.checkNoCommon().checkHaveCasesControls();
			
			
			try(VCFIterator iter1 = super.openVCFIterator(input)) {
				final VCFHeader header = iter1.getHeader();
				final OrderChecker<VariantContext> check = new OrderChecker<>(header.getSequenceDictionary(),false);

				this.casesControls.retain(header);
				this.casesControls.checkNoCommon().checkHaveCasesControls();

				
				
				final FilteringVariantContextIterator iter2 = new FilteringVariantContextIterator(iter1,CTX->{
					if(!StringUtils.isBlank(limitContig) && !limitContig.equals(CTX.getContig())) return false;
					if(CTX.getGenotypes().stream().filter(G->casesControls.contains(G)).noneMatch( G->G.hasAltAllele())) return false;
					if(!variantFilter.test(CTX))return false;
					return true;
					});
				final PeekIterator<VariantContext> iter = new PeekIterator<>(iter2);
				
				
				final SAMSequenceDictionary						
 dict = header.getSequenceDictionary();
				if(dict!=null && !StringUtils.isBlank(limitContig) && dict.getSequence(this.limitContig)==null) {
					throw new JvarkitException.ContigNotFoundInDictionary(this.limitContig, dict);
					}
				
				// create empty VCF so there is always a VCF in output
				if(this.outputVcfPath!=null ) {
					final VCFHeader hdr = new VCFHeader(header);
					try( VariantContextWriter vw=  new VariantContextWriterBuilder().
							setOutputPath(this.outputVcfPath).
							setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).
							setOption(Options.INDEX_ON_THE_FLY).
							build()) {
						vw.writeHeader(hdr);
						}
					}
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
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
					pw.print("case_alt");
					pw.print("\t");
					pw.print("case_ref");
					pw.print("\t");
					pw.print("controls_alt");
					pw.print("\t");
					pw.print("controls_ref");
					pw.print("\t");
					pw.print("variants.count");
					pw.println();
				
				
					final List<VariantContext> buffer = new ArrayList<>(100_000);
					for(;;) {
						final VariantContext ctx1 = iter.hasNext()?iter.peek():null;
						//first variant in buffer ?
						if(buffer.isEmpty() && ctx1!=null) {
							buffer.add(check.apply(iter.next()));//consumme
							continue;
							}
						// dont check for contig here !
						
						if((ctx1==null || !ctx1.contigsMatch(buffer.get(0)) || CoordMath.getLength(buffer.get(0).getStart(), ctx1.getStart()) > this.window_size) && !buffer.isEmpty()) {
							int case_alt = 0;
							int case_ref = 0;
							int ctrl_alt = 0;
							int ctrl_ref = 0;
							for(String cas : this.casesControls.getCases()) {
								boolean has_alt =  buffer.stream().map(V->V.getGenotype(cas)).anyMatch(G->G.hasAltAllele());
								if(has_alt) {
									case_alt++;
									}
								else
									{
									case_ref++;
									}
								}
							
							for(String ctrl : this.casesControls.getControls()) {
								boolean has_alt = buffer.stream().map(V->V.getGenotype(ctrl)).anyMatch(G->G.hasAltAllele());
								if(has_alt) {
									ctrl_alt++;
									}
								else
									{
									ctrl_ref++;
									}
								}
							final double pvalue = FisherExactTest.compute(case_alt, case_ref,ctrl_alt,ctrl_ref).getAsDouble();
							if(pvalue <= this.fisherTreshold) {
								final String title = new SimpleInterval(buffer.get(0).getContig(),buffer.get(0).getStart(),buffer.get(buffer.size()-1).getStart()).toString();
								pw.print(buffer.get(0).getContig());
								pw.print("\t");
								pw.print(buffer.get(0).getStart()-1);
								pw.print("\t");
								pw.print(buffer.get(buffer.size()-1).getStart());
								pw.print("\t");
								pw.print(title);
								pw.print("\t");
								pw.print(CoordMath.getLength(buffer.get(0).getStart(),buffer.get(buffer.size()-1).getStart()));
								pw.print("\t");
								pw.print(pvalue);
								pw.print("\t");
								pw.print(case_alt);
								pw.print("\t");
								pw.print(case_ref);
								pw.print("\t");
								pw.print(ctrl_alt);
								pw.print("\t");
								pw.print(ctrl_ref);
								pw.print("\t");
								pw.print(buffer.size());
								pw.println();
								
								//save vcf
								if(this.outputVcfPath!=null && pvalue < best_pvalue) {
									final VCFHeader hdr = new VCFHeader(header);
									hdr.addMetaDataLine(new VCFHeaderLine("PVALUE", String.valueOf(pvalue)));
									hdr.addMetaDataLine(new VCFHeaderLine("N_VARIANTS", String.valueOf(buffer.size())));
									try( VariantContextWriter vw=  new VariantContextWriterBuilder().
											setOutputPath(this.outputVcfPath).
											setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).
											setOption(Options.INDEX_ON_THE_FLY).
											build()) {
										vw.writeHeader(hdr);
										for(VariantContext ctx:buffer) {
											vw.add(ctx);
											}
										}
									}
								
								if( pvalue < best_pvalue) {
									LOG.info("New best p-value "+pvalue+" "+title);
									best_pvalue=pvalue;
									pw.flush();
									}
								}
	
							//reduce
							int buffer_start = buffer.get(0).getStart();
							while(!buffer.isEmpty() && CoordMath.getLength(buffer_start, buffer.get(0).getStart()) < this.window_shift) {
								buffer.remove(0);
								}
							}
						//EOF
						if(buffer.isEmpty() && ctx1==null) {
							break;
							}
					
						//consumme
						if(ctx1!=null && (buffer.isEmpty() || ctx1.contigsMatch(buffer.get(0)))) {
							buffer.add(check.apply(iter.next()));
							}
						}//end while iter
					pw.flush();
					} //end PrintWriter
				iter2.close();
				} //end VCF iterator
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

	
public static void main(final String[] args) {
	new VcfBurdenSlidingWindow().instanceMainWithExit(args);
	}
}
