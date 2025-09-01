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
package com.github.lindenb.jvarkit.tools.roh;


import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
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
	private static final Logger LOG = Logger.of(VcfROH.class);

	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;
	@Parameter(names={"--score"},description="HOM_REF/HOM_VAR score. HET score will be 'x' - 1.0.")
	private double hom_score = 0.0025;
	@Parameter(names={"--nc2hr"},description="treat NO_CALL to HOM_REF")
	private boolean nocall_to_homref = false;
	@Parameter(names={"--bed"},description="Output BED instead of interval list")
	private boolean output_as_bed = false;
	@Parameter(names={"--min-variants"},description="Minimum number of variants.")
	private int min_variants_count = 10;
	@Parameter(names={"--min-length"},description="Minimum block length. " + DistanceParser.OPT_DESCRIPTION ,converter = DistanceParser.StringConverter.class,splitter = NoSplitter.class)
	private int min_block_length = 1_000;
	@Parameter(names={"--min-score"},description="Minimum block score.")
	private double min_block_score = 0.0;
	@Parameter(names={"--regions-file"},description="Limit to intervals overlapping that BED file.")
	private Path intervalBedPath = null;


	private class Sample {
		final String name;
		final int index;
		String prev_contig = null;
		int start = -1;
		int end  = -1;
		double max_score=0;
		double score=0;
		int count = 0;
		Counter<GenotypeType> countTypes = new Counter<>();
		Sample(final String sn,int index) {
			this.name= sn;
			this.index = index;
			}
		private void dump(final PrintWriter w,final IntervalTreeMap<?> treemap) {
			if(StringUtils.isBlank(this.prev_contig)) return;
			if(this.start>=this.end) return;
			
			if(this.max_score < min_block_score) return;
			if(this.count < min_variants_count) return;
			final SimpleInterval r = new SimpleInterval(this.prev_contig, this.start,this.end);
			if( r.getLengthOnReference() < min_block_length) return;
			if(treemap!=null && !treemap.containsOverlapping(r)) {
				return;
				}
			
			w.print(r.getContig());
			w.print("\t");
			w.print(r.getStart() + (output_as_bed?-1:0));
			w.print("\t");
			w.print(r.getEnd());
			w.print("\t");
			w.print(r.getLengthOnReference());
			w.print("\t");
			w.print(this.name);
			w.print("\t");
			w.print(this.count);
			w.print("\t");
			w.print(Arrays.stream(GenotypeType.values()).filter(GT->countTypes.count(GT)>0L).map(GT->GT.name()+":"+countTypes.count(GT)).collect(Collectors.joining(";")));
			w.print("\t");
			w.print(this.max_score);
			w.println();
		}
		
		void visit(final PrintWriter w,final VariantContext ctx,final IntervalTreeMap<Boolean> treemap) {
			final Genotype gt = ctx.getGenotype(this.index);
			final boolean last_is_hom;
			// haploid
			if(gt.getAlleles().size()==1 || (gt.getAlleles().size()==2 && gt.getAlleles().stream().anyMatch(A->A.equals(Allele.SPAN_DEL)))) {
				last_is_hom = true;
				}
			else {
				switch(gt.getType()) {
					case HOM_VAR: last_is_hom =  true; break;
					case HOM_REF: last_is_hom = true; break;
					case NO_CALL: {
						last_is_hom = VcfROH.this.nocall_to_homref;
						break;
						}
					default:
						last_is_hom=false;
						break;
					}
				}
			
			final double dx;
			if(last_is_hom) {
					dx = VcfROH.this.hom_score;
				} else {
					dx = VcfROH.this.hom_score - 1.0;
				}
			
			if(!ctx.getContig().equals(this.prev_contig) || this.score+dx <= 0) {
				dump(w,treemap);
				this.prev_contig = ctx.getContig();
				this.start = ctx.getStart();
				this.score  = 0;
				this.max_score = 0;
				this.count  = 0;
				this.countTypes = new Counter<>();
				}
			this.end = ctx.getEnd();
			this.score += dx;
			this.max_score = Math.max(this.score,this.max_score);
			this.count++;
			this.countTypes.incr(gt.getType());
			}
		void finish(final PrintWriter w,final IntervalTreeMap<Boolean> treemap) {
			dump(w,treemap);
		}
	}
	
	public VcfROH() {
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			try(VCFIterator iter = input==null?new BcfIteratorBuilder().open(stdin()):new BcfIteratorBuilder().open(input)) {
				final VCFHeader header= iter.getHeader();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				final OrderChecker<VariantContext> checker = new OrderChecker<VariantContext>(dict,false);
				final List<Sample> samples = new ArrayList<>(header.getNGenotypeSamples());
				for(String sn:header.getSampleNamesInOrder()) {
					final Sample sample =new Sample(sn, header.getSampleNameToOffset().get(sn));
					samples.add(sample);
				}
				
				final IntervalTreeMap<Boolean> treemap;
				if (this.intervalBedPath!=null) {
					try(BedLineReader br = new BedLineReader(this.intervalBedPath)) {
						treemap = br.toIntervalTreeMap(B->Boolean.TRUE);
					}
				} else {
					treemap = null;
				}
				
				try(PrintWriter w =  super.openPathOrStdoutAsPrintWriter(this.output)) {
					if(!this.output_as_bed) {
						final SAMFileHeader samHeader = new SAMFileHeader(dict);
						samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
						final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
						for(Sample sn:samples) {
							final SAMReadGroupRecord rg = new SAMReadGroupRecord(sn.name);
							rg.setSample(sn.name);
							samHeader.addReadGroup(rg);
							}
						JVarkitVersion.getInstance().addMetaData(this, samHeader);
						codec.encode(w, samHeader);
						}
					while(iter.hasNext()) {
						final VariantContext ctx = checker.apply(iter.next());
						for(Sample sample:samples) sample.visit(w,ctx,treemap);
						}
					for(Sample sample:samples) sample.finish(w,treemap);
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
