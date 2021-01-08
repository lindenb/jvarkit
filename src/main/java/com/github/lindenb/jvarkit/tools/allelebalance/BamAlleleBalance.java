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
package com.github.lindenb.jvarkit.tools.allelebalance;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfDoubles;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;
/*
BEGIN_DOC

# Example

```bash
$ java -jar dist/bamallelebalance.jar --vcf src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S*.bam
RANGE	S5	S3	S2	S1	S4
[-Inf / 0.1[	0	1	1	1	1
[0.1 / 0.2[	0	0	0	0	0
[0.2 / 0.3[	0	0	0	0	0
[0.3 / 0.4[	0	0	0	0	0
[0.4 / 0.6[	0	0	0	1	1
[0.6 / 0.7[	0	0	0	0	1
[0.7 / 0.8[	0	0	0	0	2
[0.8 / 0.9[	0	0	0	0	0
[0.9 / Inf[	0	0	0	0	0
```

plotting the output

```R
T<-read.table("jeter.tsv", header=TRUE, sep="\t")
categories<-T[,1]
values<-as.matrix(T[,2:ncol(T)])
barplot(values, col = terrain.colors(length(categories)),legend=TRUE,xlab="Samples",ylab="Count",main="AD per Sample")
legend("topleft", 
       legend = categories, 
        bty = "n",
       fill = terrain.colors(length(categories))
    )
```

END_DOC
 */
@Program(name="bamallelebalance",
	description="Compute statistics about allele balance from a set of Bams",
	keywords= {"vcf","allele-balance","depth","bam"},
	modificationDate="20200805",
	creationDate="20200805"
	)
public class BamAlleleBalance extends Launcher {
	private static final Logger LOG = Logger.build(BamAlleleBalance.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx=null;
	@Parameter(names={"-v","--vcf","--variants"},description="Positions to consider. Must be frequence. Only unfiltered diallelic snp are considered.",required=true)
	private Path variantsPath = null;
	@Parameter(names={"-r","--range","--ratios"},description="Alleles Ratios. " + RangeOfDoubles.OPT_DESC ,converter=RangeOfDoubles.StringConverter.class,splitter=NoSplitter.class)
	private RangeOfDoubles ratios = new RangeOfDoubles("0.1;0.2;0.3;0.4;0.6;0.7;0.8;0.9");
	@Parameter(names={"--groupby","--partition"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"--mapq"},description="Min MAPQ")
	private int mapq = 1;
	@Parameter(names={"--min-depth"},description="Ignore sites having a depth lower than 'x'")
	private int min_depth = 10;

	private static class SiteInfo {
		int countRef = 0;
		int countAlt = 0;
	}

	private static class Snp {
		final String contigSrc;
		final int pos1;
		final char ref;
		final char alt;
		
		Snp(final VariantContext ctx) {
			this.contigSrc = ctx.getContig();
			this.pos1 = ctx.getStart();
			this.ref = ctx.getAlleles().get(0).getDisplayString().toUpperCase().charAt(0);
			this.alt = ctx.getAlleles().get(1).getDisplayString().toUpperCase().charAt(0);
		}
		
	}
	
	
	@Override
	public int doWork(List<String> args) {
		try {
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			if(bamPaths.isEmpty()) {
				LOG.error("Bam list is empty");
				return -1;
				}
			final List<Snp> snps = new ArrayList<>(10_000);
			try(VCFReader vr = VCFReaderFactory.makeDefault().open(this.variantsPath,false)) {
				try(CloseableIterator<VariantContext> iter = vr.iterator()) {
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						if(!ctx.isBiallelic() || ctx.isFiltered() || !ctx.isSNP()) continue;
						snps.add(new Snp(ctx));
						}
					}
				}
			if(snps.isEmpty()) {
				LOG.error("snp list is empty");
				return -1;
				}
			LOG.info(String.valueOf(snps.size())+" snp(s)");
			final Map<String, Counter<RangeOfDoubles.Range>> sample2count = new HashMap<>(bamPaths.size());
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			for(final Path bamPath:bamPaths) {
				try(SamReader samReader = srf.open(bamPath)) {
					final String defaultPartition= IOUtils.getFilenameWithoutCommonSuffixes(bamPath);
					final SAMFileHeader header = samReader.getFileHeader();
					final Set<String> samples = header.getReadGroups().
						stream().
						map(RG->this.groupBy.apply(RG, defaultPartition)).
						collect(Collectors.toSet())
						;
					samples.stream().forEach(SN-> sample2count.putIfAbsent(SN, new Counter<>()));
					
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
					for(final Snp snp:snps) {
						final String bamContig = converter.apply(snp.contigSrc);
						
						if(StringUtils.isBlank(bamContig)) continue;
						
						final Map<String,SiteInfo> sample2siteinfo = new HashMap<>(samples.size());
						for(final String sn :samples) {
							sample2siteinfo.put(sn, new SiteInfo());
							}
						
						try(CloseableIterator<SAMRecord> iter = samReader.queryOverlapping(bamContig,snp.pos1,snp.pos1)) {
							while(iter.hasNext()) {
								final SAMRecord rec = iter.next();
								if(!SAMRecordDefaultFilter.accept(rec, this.mapq)) continue;
								final SAMReadGroupRecord rg = rec.getReadGroup();
								if(rg==null) continue;
								String sample = this.groupBy.apply(rg, defaultPartition);
								if(StringUtils.isBlank(sample) || !sample2siteinfo.containsKey(sample)) continue;

								final Cigar cigar = rec.getCigar();
								if(cigar==null || cigar.isEmpty()) continue;
								byte bases[] = rec.getReadBases();
								
								if(bases==null || bases==SAMRecord.NULL_SEQUENCE) continue;
								
								final int readpos1 = rec.getReadPositionAtReferencePosition(snp.pos1);
								if(readpos1<1) continue;
								final int readpos0 = readpos1-1;
								if(readpos0<0 || readpos0>=bases.length) continue;
								final char base = (char)Character.toUpperCase(bases[readpos0]);
								final SiteInfo si = sample2siteinfo.get(sample);
								if(si==null) continue;
								if(base==snp.ref) {
									si.countRef++;
									}
								else if(base==snp.alt) {
									si.countAlt++;
									}
									
								}
							
						for(final String sample: sample2siteinfo.keySet()) {
							final SiteInfo si = sample2siteinfo.get(sample);
							final int depth = si.countAlt + si.countRef;
							if(depth<this.min_depth) continue;
							if(si.countAlt==0 || si.countRef==0) continue;//must be het
							final double ratio = si.countAlt/(double)depth;
							final RangeOfDoubles.Range range = this.ratios.getRange(ratio);
							sample2count.get(sample).incr(range);
							}
							
						}//end of iter
					}//end of loop over snps
				}
			}// end of loop over each bams
		
			
			// print report
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				final List<String> sortedSamples = sample2count.keySet().stream().
						sorted((A,B)->Long.compare(sample2count.get(A).getTotal(), sample2count.get(B).getTotal())).
						collect(Collectors.toList());
				
				pw.print("RANGE");
				for(final String sn: sortedSamples) {
					pw.print("\t");
					pw.print(sn);
					}
				pw.println();
				for(final RangeOfDoubles.Range r: this.ratios.getRanges()) {
					pw.print(r.toString());
					for(final String sn: sortedSamples) {
						pw.print("\t");
						pw.print(sample2count.get(sn).count(r));
						}
					pw.println();
					}
				pw.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	
	public static void main(final String[] args) {
		new BamAlleleBalance().instanceMainWithExit(args);
	}

}
