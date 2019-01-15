/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Motivation

finds the regions having some clipped reads.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
$ java -jar dist/samscansplitreads.jar src/test/resources/S*.bam 2> /dev/null 
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=M3,Number=1,Type=Float,Description="Median size of the clip in 3'">
##FORMAT=<ID=M5,Number=1,Type=Float,Description="Median size of the clip in 5'">
##FORMAT=<ID=N3,Number=1,Type=Integer,Description="Number of clipped reads in 3'">
##FORMAT=<ID=N5,Number=1,Type=Integer,Description="Number of clipped reads in 5'">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##contig=<ID=RF01,length=3302>
##contig=<ID=RF02,length=2687>
##contig=<ID=RF03,length=2592>
##contig=<ID=RF04,length=2362>
##contig=<ID=RF05,length=1579>
##contig=<ID=RF06,length=1356>
##contig=<ID=RF07,length=1074>
##contig=<ID=RF08,length=1059>
##contig=<ID=RF09,length=1062>
##contig=<ID=RF10,length=751>
##contig=<ID=RF11,length=666>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S3	S4	S5	S1	S2
RF01	195	.	N	<SPLIT>	.	.	DP=1;END=199;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:4.00:0:1	./.	./.
RF01	509	.	N	<SPLIT>	.	.	DP=1;END=577;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	725	.	N	<SPLIT>	.	.	DP=2;END=793;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF01	903	.	N	<SPLIT>	.	.	DP=2;END=1000;SVLEN=98	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	./.	./.
RF01	1607	.	N	<SPLIT>	.	.	DP=2;END=1616;SVLEN=10	GT:DP:M3:M5:N3:N5	0/1:1:0:9.00:0:1	./.	./.	./.	0/1:1:0:9.00:0:1
RF01	1672	.	N	<SPLIT>	.	.	DP=2;END=1682;SVLEN=11	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	./.	0/1:1:0:3.00:0:1	./.
RF01	1822	.	N	<SPLIT>	.	.	DP=1;END=1890;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	1926	.	N	<SPLIT>	.	.	DP=1;END=1994;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	2377	.	N	<SPLIT>	.	.	DP=2;END=2385;SVLEN=9	GT:DP:M3:M5:N3:N5	0/1:1:0:8.00:0:1	./.	./.	./.	0/1:1:0:8.00:0:1
RF01	2542	.	N	<SPLIT>	.	.	DP=2;END=2610;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:4.00:1:1	./.	./.
RF01	2689	.	N	<SPLIT>	.	.	DP=1;END=2691;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:2.00:0:1	./.	./.	./.
RF01	2719	.	N	<SPLIT>	.	.	DP=2;END=2787;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF01	3230	.	N	<SPLIT>	.	.	DP=2;END=3231;SVLEN=2	GT:DP:M3:M5:N3:N5	0/1:1:0:1.00:0:1	./.	./.	./.	0/1:1:0:1.00:0:1
RF02	3	.	N	<SPLIT>	.	.	DP=2;END=6;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF02	343	.	N	<SPLIT>	.	.	DP=4;END=451;SVLEN=109	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	0/1:2:68.00:0:2:0	0/1:1:0:3.00:0:1
RF02	513	.	N	<SPLIT>	.	.	DP=1;END=581;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	661	.	N	<SPLIT>	.	.	DP=1;END=729;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	818	.	N	<SPLIT>	.	.	DP=4;END=848;SVLEN=31	GT:DP:M3:M5:N3:N5	0/1:2:0:8.00:0:2	./.	./.	./.	0/1:2:0:8.00:0:2
RF02	957	.	N	<SPLIT>	.	.	DP=1;END=966;SVLEN=10	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	./.	./.	./.
RF02	1095	.	N	<SPLIT>	.	.	DP=1;END=1163;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	1707	.	N	<SPLIT>	.	.	DP=2;END=1725;SVLEN=19	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	0/1:1:0:2.00:0:1	./.	./.
RF02	1811	.	N	<SPLIT>	.	.	DP=1;END=1821;SVLEN=11	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:10.00:0:1	./.
RF02	1883	.	N	<SPLIT>	.	.	DP=2;END=1951;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF02	2220	.	N	<SPLIT>	.	.	DP=1;END=2224;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:4.00:0:1	./.
RF02	2515	.	N	<SPLIT>	.	.	DP=3;END=2663;SVLEN=149	GT:DP:M3:M5:N3:N5	./.	./.	0/1:3:68.00:4.00:2:1	./.	./.
RF03	500	.	N	<SPLIT>	.	.	DP=2;END=569;SVLEN=70	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	0/1:1:0:2.00:0:1	./.	./.
RF03	739	.	N	<SPLIT>	.	.	DP=2;END=823;SVLEN=85	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:2:68.00:2.00:1:1	./.
RF03	1072	.	N	<SPLIT>	.	.	DP=3;END=1147;SVLEN=76	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	0/1:1:68.00:0:1:0	./.	0/1:1:0:2.00:0:1
RF03	1207	.	N	<SPLIT>	.	.	DP=3;END=1350;SVLEN=144	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	0/1:1:68.00:0:1:0	./.
RF03	1729	.	N	<SPLIT>	.	.	DP=3;END=1809;SVLEN=81	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	0/1:1:0:7.00:0:1	./.	./.	0/1:1:68.00:0:1:0
RF03	1924	.	N	<SPLIT>	.	.	DP=2;END=1926;SVLEN=3	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	./.	./.	0/1:1:0:2.00:0:1
RF03	2153	.	N	<SPLIT>	.	.	DP=1;END=2221;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF04	173	.	N	<SPLIT>	.	.	DP=1;END=181;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:8.00:0:1	./.	./.	./.
RF04	579	.	N	<SPLIT>	.	.	DP=2;END=678;SVLEN=100	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	./.	./.
RF04	704	.	N	<SPLIT>	.	.	DP=2;END=707;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF04	754	.	N	<SPLIT>	.	.	DP=1;END=822;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF04	879	.	N	<SPLIT>	.	.	DP=1;END=887;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:8.00:0:1	./.
RF04	966	.	N	<SPLIT>	.	.	DP=2;END=1091;SVLEN=126	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	0/1:1:68.00:0:1:0	./.
RF04	1119	.	N	<SPLIT>	.	.	DP=1;END=1187;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF04	1378	.	N	<SPLIT>	.	.	DP=1;END=1380;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:2.00:0:1	./.	./.
RF04	1793	.	N	<SPLIT>	.	.	DP=2;END=1920;SVLEN=128	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:2:68.00:0:2:0	./.
RF04	2070	.	N	<SPLIT>	.	.	DP=1;END=2071;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:1.00:0:1	./.
RF05	112	.	N	<SPLIT>	.	.	DP=2;END=115;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF05	252	.	N	<SPLIT>	.	.	DP=1;END=256;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:4.00:0:1	./.
RF05	427	.	N	<SPLIT>	.	.	DP=1;END=433;SVLEN=7	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:6.00:0:1	./.	./.
RF05	529	.	N	<SPLIT>	.	.	DP=1;END=530;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:1.00:0:1	./.	./.
RF05	560	.	N	<SPLIT>	.	.	DP=2;END=563;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF05	750	.	N	<SPLIT>	.	.	DP=1;END=754;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	./.	./.
RF05	841	.	N	<SPLIT>	.	.	DP=1;END=909;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF05	960	.	N	<SPLIT>	.	.	DP=1;END=1028;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF05	1434	.	N	<SPLIT>	.	.	DP=1;END=1436;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:2.00:0:1	./.	./.	./.
RF06	26	.	N	<SPLIT>	.	.	DP=1;END=29;SVLEN=4	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:3.00:0:1	./.	./.
RF06	253	.	N	<SPLIT>	.	.	DP=1;END=321;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF06	465	.	N	<SPLIT>	.	.	DP=2;END=533;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF06	691	.	N	<SPLIT>	.	.	DP=2;END=762;SVLEN=72	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	0/1:1:68.00:0:1:0	./.
RF06	1045	.	N	<SPLIT>	.	.	DP=3;END=1133;SVLEN=89	GT:DP:M3:M5:N3:N5	./.	0/1:2:68.00:3.00:1:1	./.	0/1:1:68.00:0:1:0	./.
RF06	1224	.	N	<SPLIT>	.	.	DP=1;END=1292;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF07	223	.	N	<SPLIT>	.	.	DP=2;END=225;SVLEN=3	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	./.	./.	0/1:1:0:2.00:0:1
RF07	345	.	N	<SPLIT>	.	.	DP=2;END=420;SVLEN=76	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	0/1:1:0:4.00:0:1	./.
RF07	790	.	N	<SPLIT>	.	.	DP=1;END=792;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:2.00:0:1	./.	./.
RF07	845	.	N	<SPLIT>	.	.	DP=2;END=913;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF08	54	.	N	<SPLIT>	.	.	DP=2;END=57;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF08	295	.	N	<SPLIT>	.	.	DP=1;END=296;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:1.00:0:1	./.	./.
RF08	668	.	N	<SPLIT>	.	.	DP=1;END=736;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF08	896	.	N	<SPLIT>	.	.	DP=1;END=904;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:8.00:0:1	./.	./.	./.
RF10	192	.	N	<SPLIT>	.	.	DP=1;END=196;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	./.	./.
RF10	433	.	N	<SPLIT>	.	.	DP=1;END=501;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF11	7	.	N	<SPLIT>	.	.	DP=3;END=133;SVLEN=127	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:5.00:0:1	0/1:1:68.00:0:1:0	0/1:1:68.00:0:1:0	./.
RF11	179	.	N	<SPLIT>	.	.	DP=1;END=247;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
```
## history

 * 2018-11-14 rewritten from scratch

END_DOC
*/
@Program(name="samscansplitreads",
	description="scan split reads",
	keywords={"sam","sv","splitreads","clip"}
		)
	public class SamScanSplitReads extends Launcher {
	private static final Logger LOG = Logger.build(SamScanSplitReads.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-x","--extend"},description="extends interval by 'x' pb before merging.")
	private int extentd=10;
	@Parameter(names={"-partition","--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"-r","--rgn"},description="limit to that region CHROM:START-END")
	private String intervalStr = null;
	@Parameter(names={"-B","--bed"},description="limit to that bed file")
	private File intervalBed = null;
	@Parameter(names={"--buffer-size"},description="dump buffer every 'x' bases. Most users should not use this.")
	private int buffer_dump_size = 10_000;

	private static final byte VOID_TO_LEFT = (byte)1; 
	private static final byte RIGHT_TO_VOID = (byte)2; 
	
	private static class Arc
		{
		String sample;
		int tid;
		byte type;
		int chromStart;
		int chromEnd;
		boolean consummed = false;
		int length() {
			return chromEnd-chromStart;
			}
		}
		
	
	private void dump(
			final SAMSequenceDictionary dict,
			final IntervalTreeMap<List<Arc>> database,
			final VariantContextWriter vcw,
			final Set<String> samples,
			final Integer before
			) { 
		final Allele REF = Allele.create("N", true);
		final Allele SPLIT = Allele.create("<SPLIT>", false);
		final ContigDictComparator ctgCmp  = new ContigDictComparator(dict);
		final List<Interval> intervals  = database.keySet().stream().
				map(R-> new Interval(
					R.getContig(),
					Math.max(1,R.getStart() - this.extentd),
					R.getEnd() + this.extentd
					)).
				filter(R->(before==null?true:R.getEnd() < before.intValue())).
				sorted((A,B)->{
					int i = ctgCmp.compare(A.getContig(), B.getContig());
					if(i!=0) return i;
					i = A.getStart() - B.getStart();
					if(i!=0) return i;
					return A.getEnd() - B.getEnd();
					}).
				collect(Collectors.toList());
		
		for(final Interval interval0:intervals) {
		
			final List<Arc> arcs = database.getOverlapping(interval0).
					stream().
					flatMap(L->L.stream()).
					filter(A->!A.consummed).
					collect(Collectors.toList());
			
			if(arcs.isEmpty()) continue;
			arcs.forEach(A->A.consummed=true);
			
			
			final VariantContextBuilder vcb = new VariantContextBuilder();
			final Set<Allele> alleles = new HashSet<>();
			alleles.add(REF);
			final List<Genotype> genotypes = new ArrayList<>(samples.size());
			
			
			vcb.chr(dict.getSequence(arcs.get(0).tid).getSequenceName());
			final int chromStart = arcs.stream().mapToInt(A->A.chromStart).min().getAsInt();
			vcb.start(chromStart);
			final int chromEnd = arcs.stream().mapToInt(A->A.chromEnd).max().getAsInt();
			vcb.stop(chromEnd);
			
			vcb.attribute(VCFConstants.END_KEY, chromEnd);
			vcb.attribute("SVLEN", 1+chromEnd-chromStart);
			
			int depth = 0;
			int nsamples = 0;
			for(final String sample : samples) {
			
				final List<Arc> sampleArcs = arcs.stream().filter(A->A.sample.equals(sample)).collect(Collectors.toList());
				if(sampleArcs.isEmpty())
					{
					genotypes.add(GenotypeBuilder.createMissing(sample, 2));
					}
				else
					{
					final GenotypeBuilder gb = new GenotypeBuilder(sample);
					alleles.add(SPLIT);
					gb.alleles(Arrays.asList(REF,SPLIT));
					final int countCat1= (int)sampleArcs.stream().filter(A->A.type==VOID_TO_LEFT).count();
					final int countCat2= (int)sampleArcs.stream().filter(A->A.type==RIGHT_TO_VOID).count();
					gb.DP(countCat1+countCat2);
					gb.attribute("N5", countCat1);
					gb.attribute("N3", countCat2);
					
					if(countCat1>0)
						{
						gb.attribute("M5",Percentile.median().evaluate(sampleArcs.stream().filter(A->A.type==VOID_TO_LEFT).mapToInt(A->A.length())));
						}
					else
						{
						gb.attribute("M5",0);
						}
					if(countCat2>0) 
						{
						gb.attribute("M3", Percentile.median().evaluate(sampleArcs.stream().filter(A->A.type==RIGHT_TO_VOID).mapToInt(A->A.length())));
						}
					else
						{
						gb.attribute("M3",0);
						}
					depth+=countCat1+countCat2;
					genotypes.add(gb.make());
					++nsamples;
					}
				}
			vcb.genotypes(genotypes);
			vcb.alleles(alleles);
			vcb.attribute(VCFConstants.DEPTH_KEY, depth);
			vcb.attribute("NSAMPLES", nsamples);
			vcw.add(vcb.make());
			}
		
		}
		

		
	@Override
	public int doWork(final List<String> args) {
	 	ConcatSam.ConcatSamIterator iter = null;
	 	VariantContextWriter vcw= null;
	 	final  IntervalTreeMap<List<Arc>> database = new IntervalTreeMap<>(); 
		try {
			final ConcatSam.Factory concatFactory = new ConcatSam.Factory();
			concatFactory.setEnableUnrollList(true);
			concatFactory.addInterval(this.intervalStr);
			if(this.intervalBed!=null)
				{
				final BedLineCodec codec = new BedLineCodec();
				final BufferedReader br = IOUtils.openFileForBufferedReading(this.intervalBed);
				br.lines().
					filter(L->!StringUtil.isBlank(L)).
					map(L->codec.decode(L)).
					filter(L->L!=null).
					forEach(B->concatFactory.addInterval(B.getContig()+":"+B.getStart()+"-"+B.getEnd()));
				br.close();
				}
			
			
			iter = concatFactory.open(args);
			
			final SAMSequenceDictionary dict = iter.getFileHeader().getSequenceDictionary();
			
			final Set<String> samples = iter.getFileHeader().getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtil.isBlank(S)).collect(Collectors.toSet());
			if(samples.isEmpty())
				{
				iter.close();
				LOG.error("No samples defined");
				return -1;
				}
			
			
		
			
			final Set<VCFHeaderLine> meta=new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(meta,false,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(meta,false,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);
			meta.add(new VCFFormatHeaderLine("N5", 1, VCFHeaderLineType.Integer,"Number of clipped reads in 5'"));
			meta.add(new VCFFormatHeaderLine("N3", 1, VCFHeaderLineType.Integer,"Number of clipped reads in 3'"));
			meta.add(new VCFFormatHeaderLine("M5", 1, VCFHeaderLineType.Float,"Median size of the clip in 5'"));
			meta.add(new VCFFormatHeaderLine("M3", 1, VCFHeaderLineType.Float,"Median size of the clip in 3'"));
			meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
			meta.add(new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of sample having some split reads"));
			
			
			final VCFHeader header=new VCFHeader(meta,samples);
			JVarkitVersion.getInstance().addMetaData(this, header);
			header.setSequenceDictionary(dict);
			vcw = super.openVariantContextWriter(outputFile);
			vcw.writeHeader(header);
			
			
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.
					newInstance().
					dictionary(iter.getFileHeader().getSequenceDictionary()).
					logger(LOG).
					build();
			
			String prevContig=null;
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.apply(iter.next());
				
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				
				
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty() || !cigar.isClipped()) continue;
				
				final String sample= this.partition.getPartion(rec, null);
				if(StringUtil.isBlank(sample))continue;
				
				if(!rec.getContig().equals(prevContig)) {
					dump(dict,database,vcw,samples,null);
					database.clear();
					prevContig = rec.getContig();
					}
				else
					{
					final int before = rec.getUnclippedStart() - this.buffer_dump_size;
					dump(dict,database,vcw,samples,before);
					database.entrySet().removeIf(entries->entries.getKey().getEnd()< before);
					}

				
				
				for(int side=0;side<2;++side) {
					final int chromStart;
					final int chromEnd;
					final byte type;
					if(side==0) {
						if(!cigar.isLeftClipped())
							{
							continue;
							}
						chromStart = rec.getUnclippedStart();
						chromEnd = rec.getStart() - 1;
						type = VOID_TO_LEFT;
						}
					else  {
						if(!cigar.isRightClipped())
							{
							continue;
							}
						chromStart = rec.getStart()+1;
						chromEnd = rec.getUnclippedEnd();
						type = RIGHT_TO_VOID;
						}
					//final int length = chromEnd  - chromStart;
					//NON peux augmenter la puissance if(length < 2) continue;
					
					final Arc arc = new Arc();
					arc.sample = sample;
					arc.tid = rec.getReferenceIndex();
					arc.chromStart = chromStart;
					arc.chromEnd = chromEnd;
					arc.type = type;
					
					final Interval rgn = new Interval(rec.getReferenceName(), chromStart, chromEnd);
					List<Arc> list = database.get(rgn);
					if(list==null) {
						list = new ArrayList<>();
						database.put(rgn,list);
						}
					list.add(arc);
					}
				}
			dump(dict,database,vcw,samples,null);
			iter.close();iter=null;
			progress.close();
			vcw.close();vcw=null;
			return 0;
			} 
		catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);	
			CloserUtil.close(vcw);	
			}
		}
	
	public static void main(final String[] args) {
		new SamScanSplitReads().instanceMainWithExit(args);
	}
	// 
}
