# XContaminations

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

For @AdrienLeger2 : cross contamination between samples by looking at the homozygous genotypes.


## Usage

```
Usage: xcontaminations [options] Files
  Options:
    -factor, --factor
      Fail factor: set if (reads sample x supporting x) <= factor (reads 
      sample x supporting y)
      Default: 10
    -filter, --filter
      [20171201](moved to jexl). A JEXL Expression that will be used to filter 
      out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -ft, --frasction-treshold
      FractionTreshold treshold
      Default: 1.0E-5
    -gf, --genotype-filter
      A Java EXpression Language (JEXL) expressions to filter a genotye in a 
      VCF. Empty string will accept all genotypes. Expression returning a TRUE 
      will accept the genotypes. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -ov, --output-vcf
      output results as a vcf file; only is --sample option is set.
      Default: false
    -sample, --sample, --sample-only
      Just use sample's name. Don't use lane/flowcell/etc... data.
      Default: false
    -se, --save-every
      [20171203] In tab-delimited mode, if output file is defined save the 
      result every x seconds.
      Default: -1
    -singleton, --singleton
      [20171212] R. Redon's idea: we're not sure that the contamination comes 
      from the watched pair.. With this option, we're sure that there is only 
      one HOM_VAR on the line and no HET.
      Default: false
    -vf, --variant-filter
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * vcf
 * contamination


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew xcontaminations
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminations.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminations.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminationsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminationsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **xcontaminations** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


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


