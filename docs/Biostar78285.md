# Biostar78285

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract BAMs coverage as a VCF file.


## Usage

```
Usage: biostar78285 [options] Files
  Options:
    -B, --bed, --capture
      Limit analysis to this bed file
    -f, --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -gcw, --gc-percent-window, --gcw
      GC% window size. (if REF is defined)
      Default: 20
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --min-depth
      Min depth tresholds.
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      When using display READ_GROUPS, how should we partition the ReadGroup ? 
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -R, --reference
      Optional. Indexed fasta Reference file. This file must be indexed with 
      samtools faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * depth
 * coverage



## See also in Biostars

 * [https://www.biostars.org/p/78285](https://www.biostars.org/p/78285)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar78285
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar78285Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar78285** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ java -jar dist/biostar78285.jar -m 300   -R  ref.fa S*.bam 

##fileformat=VCFv4.2
##Biostar78285.SamFilter=record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
##FILTER=<ID=DP_LT_300,Description="All  genotypes have DP< 300">
##FORMAT=<ID=DF,Number=1,Type=Integer,Description="Number of Reads on plus strand">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of Reads on minus strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AVG_DP,Number=1,Type=Float,Description="Mean depth">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=FRACT_DP_LT_300,Number=1,Type=Float,Description="Fraction of genotypes having DP< 300">
##INFO=<ID=GC_PERCENT,Number=1,Type=Integer,Description="GC% window_size:20">
##INFO=<ID=MAX_DP,Number=1,Type=Integer,Description="Max depth">
##INFO=<ID=MEDIAN_DP,Number=1,Type=Float,Description="Median depth">
##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">
##INFO=<ID=NUM_DP_LT_300,Number=1,Type=Integer,Description="Number of genotypes having DP< 300">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	1	.	G	.	.	DP_LT_300	AVG_DP=4.25;DP=17;FRACT_DP_LT_300=1.0;GC_PERCENT=38;MAX_DP=5;MEDIAN_DP=4.50;MIN_DP=3;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:5:5:0	./.:5:5:0	./.:3:3:0	./.:4:4:0
rotavirus	2	.	G	.	.	DP_LT_300	AVG_DP=9.50;DP=38;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=14;MEDIAN_DP=8.50;MIN_DP=7;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:14:14:0	./.:9:9:0	./.:8:8:0	./.:7:7:0
rotavirus	3	.	C	.	.	DP_LT_300	AVG_DP=12.25;DP=49;FRACT_DP_LT_300=1.0;GC_PERCENT=39;MAX_DP=18;MEDIAN_DP=11.50;MIN_DP=8;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:18:18:0	./.:11:11:0	./.:12:12:0	./.:8:8:0
rotavirus	4	.	T	.	.	DP_LT_300	AVG_DP=16.25;DP=65;FRACT_DP_LT_300=1.0;GC_PERCENT=37;MAX_DP=22;MEDIAN_DP=17.00;MIN_DP=9;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:22:22:0	./.:16:16:0	./.:18:18:0	./.:9:9:0
rotavirus	5	.	T	.	.	DP_LT_300	AVG_DP=20.75;DP=83;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=27;MEDIAN_DP=21.50;MIN_DP=13;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:27:27:0	./.:18:18:0	./.:25:25:0	./.:13:13:0
rotavirus	6	.	T	.	.	DP_LT_300	AVG_DP=24.25;DP=97;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=33;MEDIAN_DP=25.50;MIN_DP=13;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:30:30:0	./.:21:21:0	./.:33:33:0	./.:13:13:0
rotavirus	7	.	T	.	.	DP_LT_300	AVG_DP=28.00;DP=112;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=38;MEDIAN_DP=30.00;MIN_DP=14;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:38:38:0	./.:23:23:0	./.:37:37:0	./.:14:14:0
rotavirus	8	.	A	.	.	DP_LT_300	AVG_DP=30.50;DP=122;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=41;MEDIAN_DP=33.00;MIN_DP=15;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:41:41:0	./.:26:26:0	./.:40:40:0	./.:15:15:0
rotavirus	9	.	A	.	.	DP_LT_300	AVG_DP=34.75;DP=139;FRACT_DP_LT_300=1.0;GC_PERCENT=44;MAX_DP=48;MEDIAN_DP=37.50;MIN_DP=16;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:46:46:0	./.:29:29:0	./.:48:48:0	./.:16:16:0
rotavirus	10	.	T	.	.	DP_LT_300	AVG_DP=40.75;DP=163;FRACT_DP_LT_300=1.0;GC_PERCENT=43;MAX_DP=56;MEDIAN_DP=43.00;MIN_DP=21;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:56:56:0	./.:35:35:0	./.:51:51:0	./.:21:21:0
rotavirus	11	.	G	.	.	DP_LT_300	AVG_DP=44.75;DP=179;FRACT_DP_LT_300=1.0;GC_PERCENT=45;MAX_DP=58;MEDIAN_DP=49.50;MIN_DP=22;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:58:58:0	./.:42:42:0	./.:57:57:0	./.:22:22:0
rotavirus	12	.	C	.	.	DP_LT_300	AVG_DP=48.75;DP=195;FRACT_DP_LT_300=1.0;GC_PERCENT=43;MAX_DP=66;MEDIAN_DP=53.00;MIN_DP=23;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:66:66:0	./.:46:46:0	./.:60:60:0	./.:23:23:0
rotavirus	13	.	T	.	.	DP_LT_300	AVG_DP=53.50;DP=214;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=73;MEDIAN_DP=58.50;MIN_DP=24;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:73:73:0	./.:51:51:0	./.:66:66:0	./.:24:24:0
```

## History

* 20180227: moved output to VCF, printing everything, adding optional BED

