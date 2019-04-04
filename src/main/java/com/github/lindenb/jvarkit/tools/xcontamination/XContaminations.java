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
package com.github.lindenb.jvarkit.tools.xcontamination;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BiPredicate;
import java.util.function.DoublePredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;
import com.github.lindenb.jvarkit.util.illumina.ShortReadName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.JexlGenotypePredicate;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

* loop over the variants in a vcf file.
* for each variant we look at the HOM (HOM_REF or HOM_VAR) variants and we look at the BAM file to test how many reads from one sample could contain the reads from another sample.

## History

* 20171122: re-written, adding support to vcf output, genotypes and variant filters.

## Input

First parameter is a VCF file or '-' for stdin.

Other parameters are a list of bam file or a file ending with '.list' and containing the path to the bam files.


## Example

```bash
$ find . -type f -name "*.bam" > bam.list
$  head -n 10000 variant.vcf | java -jar dist/xcontaminations.jar - bam.list > out.tsv
$ verticalize out.tsv


>>> 2
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ10:C3FBPACXX:0:4
$2                            sample1 : B00G5V9
$3        Machine:FlowCell:Run:Lane-2 : HISEQ10:C486PACXX:0:3
$4                            sample2 : B00G7LK
$5                          same.lane : 0
$6   reads_sample1_supporting_sample1 : 26392
$7   reads_sample1_supporting_sample2 : 70
$8    reads_sample1_supporting_others : 40
$9   reads_sample2_supporting_sample2 : 21473
$10  reads_sample2_supporting_sample1 : 39
$11    reads_sample2_supporting_other : 31
<<< 2

(...)

>>> 9
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ5:C3FV0ACXX:0:7
$2                            sample1 : B00G738
$3        Machine:FlowCell:Run:Lane-2 : HISEQ5:C3FV0ACXX:0:7
$4                            sample2 : B00G754
$5                          same.lane : 1
$6   reads_sample1_supporting_sample1 : 10209
$7   reads_sample1_supporting_sample2 : 23
$8    reads_sample1_supporting_others : 15
$9   reads_sample2_supporting_sample2 : 9054
$10  reads_sample2_supporting_sample1 : 32
$11    reads_sample2_supporting_other : 9
<<< 9

```

generating in parallel:

```make
CHROMS=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

.PHONY:all

define xcont

$$(addprefix tmp.xcont.,$$(addsuffix .tsv.gz,$(1))) :
        bcftools view -r "$(1)" unifiedgenotyper.vcf.gz -Tcapture.bed | java -Xmx1g -jar xcontaminations.jar - bal.list  | gzip --best > $$(addsuffix .tmp.gz,$$@) && mv  
$$(addsuffix .tmp.gz,$$@) $$@


endef

all: $(foreach C,${CHROMS},$(addprefix tmp.xcont.,$(addsuffix .tsv.gz,${C})))
        $(foreach I,0 1, gunzip -c $^ |  awk -F '       ' '($$5==$I)'  |awk -F '        ' 'BEGIN {T=0;N=0;} {for(i=6;i<=NF;++i) T+=int($$i); N+=int($$7); N+=int($$10);} E
ND { printf("%f\n",N/T);}'; )

$(foreach C,${CHROMS},$(eval $(call xcont,$C)))
```

## Example

vcf output:


```
$ java -jar dist/xcontaminations.jar -ov -sample -vf 'DP>100' mutations.vcf *.bam 
```



```
##fileformat=VCFv4.2
##FILTER=<ID=BADSAMPLES,Description="At least one pair of genotype fails the 'LE' test">
##FILTER=<ID=XCONTAMINATION,Description="Fisher test is < 1.0E-5">
##FORMAT=<ID=F,Number=1,Type=Float,Description="Fisher test. '-1' for unavailable.">
##FORMAT=<ID=S1A,Number=1,Type=Character,Description="sample 1 allele">
##FORMAT=<ID=S1S1,Number=1,Type=Integer,Description="reads sample 1 supporting sample 1">
##FORMAT=<ID=S1S2,Number=1,Type=Integer,Description="reads sample 1 supporting sample 2">
##FORMAT=<ID=S1SO,Number=1,Type=Integer,Description="reads sample 1 supporting others">
##FORMAT=<ID=S2A,Number=1,Type=Character,Description="sample 2 allele">
##FORMAT=<ID=S2S1,Number=1,Type=Integer,Description="reads sample 2 supporting sample 1">
##FORMAT=<ID=S2S2,Number=1,Type=Integer,Description="reads sample 2 supporting sample 2">
##FORMAT=<ID=S2SO,Number=1,Type=Integer,Description="reads sample 2 supporting others">
##INFO=<ID=BADSAMPLES,Number=.,Type=String,Description="Samples founds failing the 'LE' test">
##INFO=<ID=LE,Number=1,Type=Integer,Description="number of pair of genotypes having (S1S1<=S1S2 or S2S2<=S2S1).">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1:S2	S1:S3	S1:S4	S2:S3	S2:S4	S3:S4
rotavirus	51	.	A	G	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S4;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:A:235:0:19:G:71:0:9	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:A:204:0:14:G:71:0:9	1.00:A:261:0:13:G:71:0:9
rotavirus	536	.	A	T	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S1;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	1.00:T:20:505:29:A:21:542:20	0.880:T:20:505:29:A:26:692:20	0.531:T:20:505:29:A:10:189:8	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	693	.	T	G	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S1;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	0.884:G:26:326:1:T:25:294:1	0.892:G:26:326:1:T:33:432:4	0.528:G:26:326:1:T:6:106:1	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	799	.	A	C	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S3;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:A:273:0:24:C:420:0:31	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:A:291:0:29:C:420:0:31	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:C:0:420:31:A:0:86:8
rotavirus	812	.	G	T	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S3;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:283:0:13:T:443:0:29	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:291:0:25:T:443:0:29	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:T:0:443:29:G:0:90:3
rotavirus	833	.	G	A	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S1;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	1.00:A:0:261:24:G:0:302:26	1.00:A:0:261:24:G:0:430:25	1.00:A:0:261:24:G:0:85:4	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	916	.	A	T	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S1,S4;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	0.188:T:16:295:0:A:23:269:0	0.530:T:16:295:0:A:28:405:0	0.091:T:16:295:0:A:10:86:0	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	1044	.	A	T	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S2,S3;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	0.265:A:123:7:0:T:144:15:0	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	0.712:T:15:144:0:A:17:139:0	0.251:T:15:144:0:A:2:56:0	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	1045	.	C	G	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S3;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:C:118:0:8:G:145:0:9	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:C:139:0:12:G:145:0:9	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:0:145:9:C:0:54:3
rotavirus	1054	.	C	G	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S2;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	1.00:C:82:0:5:G:97:0:5	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:0:97:5:C:0:88:11	1.00:G:0:97:5:C:0:39:1	-1.0:.:-1:-1:-1:.:-1:-1:-1
rotavirus	1064	.	G	A	.	BADSAMPLES;XCONTAMINATION	BADSAMPLES=S4;LE=3	F:S1A:S1S1:S1S2:S1SO:S2A:S2S1:S2S2:S2SO	-1.0:.:-1:-1:-1:.:-1:-1:-1	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:47:0:1:A:24:0:0	-1.0:.:-1:-1:-1:.:-1:-1:-1	1.00:G:57:0:6:A:24:0:0	1.00:G:49:0:4:A:24:0:0
```


END_DOC

 */

@Program(name="xcontaminations",
	description="For @AdrienLeger2 : cross contamination between samples by looking at the homozygous genotypes.",
	keywords= {"sam","bam","vcf","contamination"}
	)
public class XContaminations extends Launcher
	{
	private static final Logger LOG=Logger.build(XContaminations.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-filter","--filter"},description="[20171201](moved to jexl). "+SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-sample","--sample","--sample-only"},description="Just use sample's name. Don't use lane/flowcell/etc... data.")
	private boolean use_only_sample_name = false;
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-gf","--genotype-filter"},description=JexlGenotypePredicate.PARAMETER_DESCRIPTION,converter=JexlGenotypePredicate.Converter.class)
	private BiPredicate<VariantContext,Genotype> genotypeFilter = JexlGenotypePredicate.create("");
	@Parameter(names={"-ov","--output-vcf"},description="output results as a vcf file; only is --sample option is set.")
	private boolean output_as_vcf= false;
	@Parameter(names={"-ft","--frasction-treshold"},description="FractionTreshold treshold")
	private double fraction_treshold=1E-5;
	@Parameter(names={"-factor","--factor"},description="Fail factor: set if (reads sample x supporting x) <= factor (reads sample x supporting y)")
	private int fail_factor=10;
	@Parameter(names={"-se","--save-every"},description="[20171203] In tab-delimited mode, if output file is defined save the result every x seconds.")
	private long save_every_sec = -1;
	@Parameter(names={"-singleton","--singleton"},description="[20171212] R. Redon's idea: we're not sure that the contamination comes from the watched pair."
			+ ". With this option, we're sure that there is only one HOM_VAR on the line and no HET.")
	private boolean use_singleton = false;

	
	private DoublePredicate passFractionTreshold  = (V) -> V > fraction_treshold;
	
	private static class SampleAlleles
		{
		long reads_sample1_supporting_sample1 = 0L;
		long reads_sample1_supporting_sample2 = 0L;
		long reads_sample1_supporting_other = 0L;
		long reads_sample2_supporting_sample1 = 0L;
		long reads_sample2_supporting_sample2 = 0L;
		long reads_sample2_supporting_other = 0L;
		long number_of_comparaisons = 0L;
		
		public double getFraction() {
			final double t= 
					  reads_sample1_supporting_sample1 +
					  reads_sample1_supporting_sample2 +
					  reads_sample2_supporting_sample1 +
					  reads_sample2_supporting_sample2
					  ;
			final double n =  
					reads_sample1_supporting_sample2 +
					reads_sample2_supporting_sample1 
					;
			return n / t;
			}
		
		@SuppressWarnings("unused")
		public double getFisher() {
			final FisherExactTest fisher = FisherExactTest.compute(
					(int)this.reads_sample1_supporting_sample1,
					(int)this.reads_sample1_supporting_sample2,
					(int)this.reads_sample2_supporting_sample1,
					(int)this.reads_sample2_supporting_sample2
					);
			return fisher.getAsDouble();
			}
		}
	
	private static class SamplePair
		{
		final SampleIdentifier sample1;
		final SampleIdentifier sample2;
		final int _hash;
		SamplePair(
				final SampleIdentifier s1,
				final SampleIdentifier s2
				)
			{
			if(s1.getSampleName().compareTo(s2.getSampleName())<0)
				{
				this.sample1=s1;
				this.sample2=s2;
				}
			else
				{
				this.sample1=s2;
				this.sample2=s1;
				}
			final int prime = 31;
			int result = 1;
			result = prime * result + sample1.hashCode();
			result = prime * result + sample2.hashCode();
			this._hash = result;
			}
		@Override
		public int hashCode() {
			return this._hash;
			}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			final SamplePair other = (SamplePair) obj;
			if(this._hash!=other._hash) return false;
			if (!sample1.equals(other.sample1)) return false;
			if (!sample2.equals(other.sample2)) return false;
			return true;
			}
		/** used as sample's name in VCF header for VCF output */
		public String getLabel() {
			return sample1.getSampleName()+":"+sample2.getSampleName();
		}
		
		@Override
		public String toString() {
			return "("+sample1+"/"+sample2+")";
			}
		}
	
	private static interface SampleIdentifier
		{
		public String getSampleName();
		public String getLabel();
		}
	
	private static class SimpleSampleIdenfifier
		implements SampleIdentifier
		{
		private final String sampleName;
		SimpleSampleIdenfifier(final String sampleName) {
			this.sampleName= sampleName;
			}
		@Override
		public String getLabel() {
			return sampleName;
			}
		@Override
		public String getSampleName() {
			return sampleName;
			}
		@Override
		public int hashCode() {
			return sampleName.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			return sampleName.equals(SimpleSampleIdenfifier.class.cast(obj).sampleName);
			}
		@Override
		public String toString() {
			return sampleName;
			}
		}
	
	private static class SequencerFlowCellRunLaneSample
		implements SampleIdentifier
		{
		private final String machine;
		private final String flowCell;
		private final int run;
		private final int lane;
		private final String sampleName;
		private final int _hash;
		SequencerFlowCellRunLaneSample(final ShortReadName name,final String sampleName)
			{
			this.machine=name.getInstrumentName();
			this.flowCell=name.getFlowCellId();
			this.run=Math.max(name.getRunId(),0);
			this.lane=name.getFlowCellLane();
			this.sampleName=sampleName;
			
			final int prime = 31;
			int result = 1;
			result = prime * result + flowCell.hashCode();
			result = prime * result + lane;
			result = prime * result + machine.hashCode();
			result = prime * result + sampleName.hashCode();
			result = prime * result + run;
			this._hash = result;
			}

		@Override
		public int hashCode() {
			return  _hash;
			}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			final SequencerFlowCellRunLaneSample other = (SequencerFlowCellRunLaneSample) obj;
			if (lane != other.lane) return false;
			if (run != other.run) return false;
			if (!flowCell.equals(other.flowCell)) return false;
			if (!machine.equals(other.machine)) return false;
			if (!sampleName.equals(other.sampleName)) return false;
			return true;
			}
		
		@Override
		public String getSampleName() {
			return this.sampleName;
			}
		
		@Override
		public String getLabel()
			{
			return machine+":"+flowCell+":"+run+":"+lane;
			}
		
		@Override
		public String toString() {
			return getLabel()+":"+sampleName;
			}
		}
	
	
	private void saveToFile(final Map<SamplePair,SampleAlleles> contaminationTable) throws IOException{
		PrintWriter pw = null;
		try 
			{
			boolean somethingPrinted=false;
			pw= super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			/* we're done, print the result */
			pw.print("#");
			if(!this.use_only_sample_name) {
				pw.print("Machine:FlowCell:Run:Lane-1");
				pw.print('\t');
				}
			pw.print("sample1");
			pw.print('\t');
			if(!this.use_only_sample_name) {
				pw.print("Machine:FlowCell:Run:Lane-2");
				pw.print('\t');
				}
			pw.print("sample2");
			pw.print('\t');
			if(!this.use_only_sample_name) {
				pw.print("same.lane");
				pw.print('\t');
				}
			pw.print("reads_sample1_supporting_sample1");
			pw.print('\t');
			pw.print("reads_sample1_supporting_sample2");
			pw.print('\t');
			pw.print("reads_sample1_supporting_others");
			pw.print('\t');
			pw.print("reads_sample2_supporting_sample2");
			pw.print('\t');
			pw.print("reads_sample2_supporting_sample1");
			pw.print('\t');
			pw.print("reads_sample2_supporting_other");
			pw.print('\t');
			pw.print("Fraction");
			pw.print('\t');
			pw.print("Pass-Fraction");
			pw.print('\t');
			pw.print("count_comparison");
			pw.println();
			for(final SamplePair pair : contaminationTable.keySet())
				{
				final SampleAlleles sampleAlleles = contaminationTable.get(pair);
				if(sampleAlleles==null) continue;
				
				if(!this.use_only_sample_name) {
					pw.print(pair.sample1.getLabel());
					pw.print('\t');
					}
				pw.print(pair.sample1.getSampleName());
				pw.print('\t');
				if(!this.use_only_sample_name) {
					pw.print(pair.sample2.getLabel());
					pw.print('\t');
					}
				pw.print(pair.sample2.getSampleName());
				pw.print('\t');
				if(!this.use_only_sample_name) {
					pw.print(pair.sample1.getLabel().equals(pair.sample2.getLabel())?1:0);
					pw.print('\t');
					}
				pw.print(sampleAlleles.reads_sample1_supporting_sample1);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample1_supporting_sample2);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample1_supporting_other);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_sample2);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_sample1);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_other);
				pw.print('\t');
				final double fraction = sampleAlleles.getFraction();
				pw.print(fraction);
				pw.print('\t');
				pw.print(this.passFractionTreshold.test(fraction)?".":"*");
				pw.print('\t');
				pw.print(sampleAlleles.number_of_comparaisons);
				pw.println();
				somethingPrinted=true;				
				}
			pw.flush();
			pw.close();
			pw=null;
			
			if(!somethingPrinted)
				{
				LOG.warn("Warning: NO output");
				}
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		long last_save_ms = System.currentTimeMillis();
		if(this.output_as_vcf && !this.use_only_sample_name)
			{
			LOG.error("cannot write vcf if --sample is not set");
			return -1;
			}
		if(args.size()<2)
			{
			LOG.error("Illegal Number of args");
			return -1;
			}
		final Set<File> bamFiles = IOUtils.unrollFiles(args.subList(1, args.size())).
				stream().map(S->new File(S)).collect(Collectors.toSet());
			
		if(bamFiles.isEmpty())
			{
			LOG.error("Undefined BAM file(s)");
			return -1;
			}	
		
		SAMRecordIterator iter=null;
		VCFIterator in=null;
		Map<String,SamReader> sample2samReader=new HashMap<>();
		VariantContextWriter vcfw = null;
		try {
			final SamReaderFactory srf= super.createSamReaderFactory();
			
			if(args.get(0).equals("-"))
				{
				in = super.openVCFIterator(null);
				}
			else
				{
				in = super.openVCFIterator(args.get(0));
				}
			
			
			VCFHeader vcfHeader=in.getHeader();
			final SAMSequenceDictionary dict1=vcfHeader.getSequenceDictionary();
			if(dict1==null)
				{
				LOG.error(JvarkitException.VcfDictionaryMissing.getMessage(args.get(0)));
				return -1;
				}
			
			final Set<String> sampleNames= new HashSet<>(vcfHeader.getSampleNamesInOrder());
			if( sampleNames.isEmpty())
				{
				LOG.error("VCF contains no sample");
				return -1;
				}
			
			for(final File bamFile:bamFiles)
				{
				LOG.info("Opening "+bamFile);
				final SamReader samReader=srf.open(bamFile);
				final SAMFileHeader samHeader= samReader.getFileHeader();
				final SAMSequenceDictionary dict2=samHeader.getSequenceDictionary();
				if(dict2==null)
					{
					samReader.close();
					LOG.error(JvarkitException.BamDictionaryMissing.getMessage(bamFile.getPath()));
					return -1;
					}
				
				if(!SequenceUtil.areSequenceDictionariesEqual(dict1, dict2))
					{
					samReader.close();
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict1, dict2));
					return -1;
					}
				
				if(!samReader.hasIndex())
					{
					samReader.close();
					LOG.error("sam is not indexed : "+bamFile);
					return -1;
					}
				String sampleName=null;
				for(final SAMReadGroupRecord rgr:samHeader.getReadGroups())
					{
					final String s=rgr.getSample();
					if(StringUtil.isBlank(s)) continue;
					if(sampleName==null)
						{
						sampleName=s;
						}
					else if(!sampleName.equals(s))
						{
						samReader.close();
						LOG.error("Cannot handle more than one sample/bam  "+bamFile+" "+sampleName);
						return -1;
						}
					}
				if(sampleName==null)
					{
					samReader.close();
					LOG.error("No sample in "+bamFile);
					continue;//skip this bam
					}
				if(!sampleNames.contains(sampleName))
					{
					samReader.close();
					LOG.error("Not in VCF header: sample "+sampleName+" "+bamFile);
					continue;//skip this bam
					}
				if(sample2samReader.containsKey(sampleName))
					{
					samReader.close();
					LOG.error("Cannot handle more than one bam/sample: "+bamFile+" "+sampleName);
					return -1;
					}
				
				sample2samReader.put(sampleName, samReader);
				}
			
			if(sample2samReader.size()<2)
				{
				LOG.error("Not engough BAM/samples. Expected at least two valid BAMs");
				return -1;
				}
			
			sampleNames.retainAll(sample2samReader.keySet());
			

			/* create a VCF is VCF output asked */
			
			final List<SamplePair> sampleListForVcf;
			if(this.output_as_vcf)
				{
				vcfw = super.openVariantContextWriter(outputFile);
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				metaData.add(new VCFFormatHeaderLine("S1S1", 1, VCFHeaderLineType.Integer,"reads sample 1 supporting sample 1"));
				metaData.add(new VCFFormatHeaderLine("S1S2", 1, VCFHeaderLineType.Integer,"reads sample 1 supporting sample 2"));
				metaData.add(new VCFFormatHeaderLine("S1SO", 1, VCFHeaderLineType.Integer,"reads sample 1 supporting others"));
				
				metaData.add(new VCFFormatHeaderLine("S2S1", 1, VCFHeaderLineType.Integer,"reads sample 2 supporting sample 1"));
				metaData.add(new VCFFormatHeaderLine("S2S2", 1, VCFHeaderLineType.Integer,"reads sample 2 supporting sample 2"));
				metaData.add(new VCFFormatHeaderLine("S2SO", 1, VCFHeaderLineType.Integer,"reads sample 2 supporting others"));
				metaData.add(new VCFFormatHeaderLine("FR", 1, VCFHeaderLineType.Float,"Fraction. '-1' for unavailable."));
				
				metaData.add(new VCFFormatHeaderLine("S1A", 1, VCFHeaderLineType.Character,"sample 1 allele"));
				metaData.add(new VCFFormatHeaderLine("S2A", 1, VCFHeaderLineType.Character,"sample 2 allele"));

				
				metaData.add(new VCFFilterHeaderLine("XCONTAMINATION","Fraction test is > "+fraction_treshold));
				metaData.add(new VCFFilterHeaderLine("BADSAMPLES","At least one pair of genotype fails the 'LE' test"));
				metaData.add(new VCFInfoHeaderLine("LE",1, VCFHeaderLineType.Integer,"number of pair of genotypes having (S1S1<=S1S2 or S2S2<=S2S1)."));
				metaData.add(new VCFInfoHeaderLine("BADSAMPLES",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"Samples founds failing the 'LE' test"));
				
				
				
				sampleListForVcf = new ArrayList<>();
				final List<String> sampleList=new ArrayList<>(sampleNames);
				for(int x=0;x+1<sampleList.size();++x) {
					for(int y=x+1;y<sampleList.size();++y) {
						sampleListForVcf.add(new SamplePair(new SimpleSampleIdenfifier(sampleList.get(x)),new SimpleSampleIdenfifier(sampleList.get(y))));
					}
				}
				
				final VCFHeader header2 = new VCFHeader(metaData, sampleListForVcf.stream().
						map(V->V.getLabel()).
						sorted().collect(Collectors.toList())
						);
				header2.setSequenceDictionary(dict1);
				vcfw.writeHeader(header2);
				}
			else
				{
				vcfw = null;
				sampleListForVcf = null;
				}
			
			final Map<SamplePair,SampleAlleles> contaminationTable=new HashMap<>();
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1).logger(LOG);
			while(in.hasNext())
				{
				final VariantContext ctx= progress.watch(in.next());
				if(!ctx.isSNP() || ctx.isFiltered() || !ctx.isBiallelic() || ctx.isSymbolic() || !this.variantFilter.test(ctx)) {
					continue;
				}
				
				int count_homref=0;
				int count_homvar=0;
				int count_het=0;
				
				final Map<String,Genotype> sample2gt = new HashMap<>();
				for(int gidx=0;gidx < ctx.getNSamples();++gidx) {
					final Genotype G = ctx.getGenotype(gidx);
					if(!G.isCalled()) continue;
					if(G.isHet())
						{
						count_het++;// here because in use_singleton we must be sure that there is only one hom_var
						if(this.use_singleton && count_het>0) break;
						}
					else if(G.isHomVar())
						{
						count_homvar++;// here because in use_singleton we must be sure that there is only one hom_var
						if(this.use_singleton && count_homvar>1) break;
						}
					
					if(G.isFiltered()) continue;
					if(!sample2samReader.containsKey(G.getSampleName())) continue;
					if(!sampleNames.contains(G.getSampleName())) continue;
					if(!this.genotypeFilter.test(ctx, G)) continue;
					sample2gt.put(G.getSampleName(), G);
				}
				if(this.use_singleton && count_het>0) continue;
				if(this.use_singleton && count_homvar>1) continue;
				
				if(sample2gt.size()<2) continue;
				
				
				//reset and recount
				count_homref =0;
				count_homvar =0;
				count_het = 0;
				for(final String sampleName:sample2gt.keySet()) {
					final Genotype G = ctx.getGenotype(sampleName);
					switch(G.getType()) {
						case HOM_REF :  count_homref++;break;
						case HOM_VAR :  count_homvar++;break;
						case HET :  count_het++;break;
						default:break;
						}
					}
				
				
								
				// singleton check
				if(this.use_singleton && ( count_het>0 || count_homvar!=1 ))
					{
					continue;
					}
				//at least one HOM_REF and one HOM_VAR
				if(count_homref==0) continue;
				if(count_homvar==0) continue;
				
						
				final Map<SampleIdentifier,Counter<Character>> sample_identifier_2allelesCount=new HashMap<>();
				
				/* scan Reads for those Genotype/Samples */
				for(final String sampleName: sample2gt.keySet())
					{
					if(!sample2samReader.containsKey(sampleName)) continue;
					//sample name is not in vcf header
					final SamReader samReader = sample2samReader.get(sampleName);
					if(samReader==null) continue;
					
					final Genotype genotype = sample2gt.get(sampleName);
					if(genotype==null) continue;
					
					iter = samReader.query(
							ctx.getContig(),
							ctx.getStart(),
							ctx.getEnd(),
							false
							);
					while(iter.hasNext())
						{
						final SAMRecord record= iter.next();
						if(record.getEnd()< ctx.getStart()) continue;
						if(ctx.getEnd()< record.getStart()) continue;
						
						if(record.getReadUnmappedFlag()) continue;
						if(this.filter.filterOut(record)) continue;
					
						
						final SAMReadGroupRecord srgr = record.getReadGroup();
						//not current sample
						if(srgr==null) continue;
						if(!sampleName.equals(srgr.getSample())) continue;
						
						final Cigar cigar=record.getCigar();
						if(cigar==null || cigar.isEmpty()) continue;
						byte readSeq[]=record.getReadBases();
						if(readSeq==null || readSeq.length==0) continue;
						
						int readPos = record.getReadPositionAtReferencePosition(ctx.getStart());
						if(readPos<1) continue;
						readPos--;
						if(readPos>=readSeq.length) continue;
						final char base = Character.toUpperCase((char)readSeq[readPos]);
						
						if(base=='N') continue;
						
						
						final SampleIdentifier sampleIdentifier;
						if(this.use_only_sample_name)
							{
							sampleIdentifier = new SimpleSampleIdenfifier(sampleName);
							}
						else
							{
							final ShortReadName readName = ShortReadName.parse(record);
							if(!readName.isValid())
								{
								LOG.info("No a valid read name "+record.getReadName());
								continue;
								}
							sampleIdentifier = new SequencerFlowCellRunLaneSample(readName, sampleName);
							}

						
						Counter<Character> sampleAlleles= sample_identifier_2allelesCount.get(sampleIdentifier);
						if(sampleAlleles==null)
							{
							sampleAlleles=new Counter<Character>();
							sample_identifier_2allelesCount.put(sampleIdentifier, sampleAlleles);
							}
						sampleAlleles.incr(base);
						}
					iter.close();
					iter=null;
					}/* end scan reads for this sample */
				
				
			
				
				/* sum-up data for this SNP */
				final VariantContextBuilder vcb;
				final List<Genotype> genotypeList;
				
				if(this.output_as_vcf)
					{
					vcb = new VariantContextBuilder(args.get(0), ctx.getContig(), ctx.getStart(), ctx.getEnd(), ctx.getAlleles());
					if(ctx.hasID()) vcb.id(ctx.getID());
					genotypeList= new ArrayList<>();
					}
				else
					{
					vcb = null;
					genotypeList = null;
					}
				
				
				for(final String sample1: sample2gt.keySet())
					{
					final Genotype g1= sample2gt.get(sample1);
					final char a1 = g1.getAllele(0).getBaseString().charAt(0);
					
					
					for(final String sample2:  sample2gt.keySet())
						{
						if(sample1.compareTo(sample2)>=0) continue;
						final Genotype g2= sample2gt.get(sample2);
						if(g2.sameGenotype(g1)) continue;
						final char a2 =  g2.getAllele(0).getBaseString().charAt(0);
						
						for(final SampleIdentifier sfcr1: sample_identifier_2allelesCount.keySet())
							{
							if(!sfcr1.getSampleName().equals(sample1)) continue;
							final Counter<Character> counter1 =  sample_identifier_2allelesCount.get(sfcr1);
							if(counter1==null) continue;

							
							for(final SampleIdentifier sfcr2: sample_identifier_2allelesCount.keySet())
								{
								if(!sfcr2.getSampleName().equals(sample2)) continue;
								
								final SamplePair samplePair = new SamplePair(sfcr1, sfcr2);
								
								final Counter<Character> counter2 =  sample_identifier_2allelesCount.get(sfcr2);
								if(counter2==null) continue;
								
								
								SampleAlleles sampleAlleles = contaminationTable.get(samplePair);
								if(sampleAlleles==null)
									{
									sampleAlleles=new SampleAlleles();
									contaminationTable.put(samplePair,sampleAlleles);
									if(!this.output_as_vcf && contaminationTable.size()%10000==0) LOG.info("n(pairs)=" + contaminationTable.size() ); 
									}

								sampleAlleles.number_of_comparaisons++;
								
								for(final Character allele: counter1.keySet())
									{
									final long n = counter1.count(allele);
									if(allele.equals(a1))
										{
										sampleAlleles.reads_sample1_supporting_sample1 += n;
										}
									else if(allele.equals(a2))
										{
										sampleAlleles.reads_sample1_supporting_sample2 += n;
										}
									else
										{
										sampleAlleles.reads_sample1_supporting_other += n;
										}
									}
								
								for(final Character allele: counter2.keySet())
									{
									final long n = counter2.count(allele);
									if(allele.equals(a2))
										{
										sampleAlleles.reads_sample2_supporting_sample2 += n;
										}
									else if(allele.equals(a1))
										{
										sampleAlleles.reads_sample2_supporting_sample1 += n;
										}
									else
										{
										sampleAlleles.reads_sample2_supporting_other += n;
										}
									}
								}
							}
						}
					}
				
				if(this.output_as_vcf) 
					{
					final Set<String> bad_samples=new TreeSet<>();
					boolean fraction_flag=false;
					int num_lt=0;
					for(final SamplePair samplepair :sampleListForVcf)
						{
						final GenotypeBuilder gb = new GenotypeBuilder(samplepair.getLabel());
						final SampleAlleles sampleAlleles = contaminationTable.get(samplepair);
						if(sampleAlleles != null)
							{
							gb.attribute("S1S1", sampleAlleles.reads_sample1_supporting_sample1);
							gb.attribute("S1S2", sampleAlleles.reads_sample1_supporting_sample2);
							gb.attribute("S1SO", sampleAlleles.reads_sample1_supporting_other);
							gb.attribute("S2S1", sampleAlleles.reads_sample2_supporting_sample1);
							gb.attribute("S2S2", sampleAlleles.reads_sample2_supporting_sample2);
							gb.attribute("S2SO", sampleAlleles.reads_sample2_supporting_other);
							gb.attribute("S1A",sample2gt.get(samplepair.sample1.getSampleName()).getAllele(0).getDisplayString().charAt(0));
							gb.attribute("S2A",sample2gt.get(samplepair.sample2.getSampleName()).getAllele(0).getDisplayString().charAt(0));
							final double fraction = sampleAlleles.getFraction();
							gb.attribute("FR", fraction);
							if(!this.passFractionTreshold.test(fraction)) {
								fraction_flag=true;
								}
							
							boolean bad_lt_flag=false;
							if( sampleAlleles.reads_sample1_supporting_sample1 <= this.fail_factor*sampleAlleles.reads_sample1_supporting_sample2)
								{
								bad_samples.add(samplepair.sample1.getSampleName());
								bad_lt_flag = true;
								}
							if(sampleAlleles.reads_sample2_supporting_sample2 <= this.fail_factor*sampleAlleles.reads_sample2_supporting_sample1) {
								bad_samples.add(samplepair.sample2.getSampleName());
								bad_lt_flag = true;
								}
							
							if(bad_lt_flag)
								{
								num_lt++;
								}
							}
						else
							{
							gb.attribute("S1S1", -1);
							gb.attribute("S1S2", -1);
							gb.attribute("S1SO", -1);
							gb.attribute("S2S1", -1);
							gb.attribute("S2S2", -1);
							gb.attribute("S2SO", -1);
							gb.attribute("S1A",'.');
							gb.attribute("S2A",'.');

							gb.attribute("FR", -1f);
							}
						genotypeList.add(gb.make());
						}
					if(!bad_samples.isEmpty())
						{
						vcb.attribute("BADSAMPLES", new ArrayList<>(bad_samples));
						}
					vcb.attribute("LE", num_lt);
					if(fraction_flag || !bad_samples.isEmpty()) 
						{
						if(fraction_flag) vcb.filter("XCONTAMINATION");
						if(!bad_samples.isEmpty()) vcb.filter("BADSAMPLES");
						}
					else
						{
						vcb.passFilters();
						}
					vcb.genotypes(genotypeList);
					vcfw.add(vcb.make());
					
					contaminationTable.clear();
					}
				else
					{
					final long now=System.currentTimeMillis();
					if(	this.outputFile!=null && 
						this.save_every_sec>-1L && 
						last_save_ms+(this.save_every_sec*1000L)> now
						) {
						saveToFile(contaminationTable);
						last_save_ms = now;
						}
					}
				}
			progress.finish();
						
			if(this.output_as_vcf)
				{
				vcfw.close();
				vcfw=null;
				}
			else
				{
				saveToFile(contaminationTable);
				}
			return 0;
			}
		catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfw);			
			CloserUtil.close(in);
			CloserUtil.close(iter);
			for(SamReader samReader:sample2samReader.values())
				CloserUtil.close(samReader);
			sample2samReader.clear();
			}
		
		}

	
	
	public static void main(final String[] args)
		{
		new XContaminations().instanceMainWithExit(args);
		}
	}
