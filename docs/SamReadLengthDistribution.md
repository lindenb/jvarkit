# SamReadLengthDistribution

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sam read/insert length distribution


## Usage

```
Usage: samreadlengthdistribution [options] Files
  Options:
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -Q, --mapq
      min MAPQ
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --version
      print version and exit
    -w, --windows
      Histogram windows. A 'range of integers' is a list of integers in 
      ascending order separated with semicolons.
      Default: [[-Inf/0[, [0/10[, [10/20[, [20/30[, [30/40[, [40/50[, [50/100[, [100/150[, [150/200[, [200/250[, [250/300[, [300/350[, [350/400[, [400/450[, [450/500[, [500/1000[, [1000/Inf[]
    -m
      method, how should I get the read/insert length ?
      Default: SEQ_LENGTH
      Possible Values: [SEQ_LENGTH, CIGAR_REF_LENGTH, CIGAR_PADDED_REF_LENGTH, INSERT_LENGTH]

```


## Keywords

 * sam
 * bam
 * histogram


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samreadlengthdistribution
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20160922

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistribution.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistribution.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistributionTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistributionTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samreadlengthdistribution** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Because gatk is buggy: http://gatkforums.broadinstitute.org/gatk/discussion/8342/duplicate-columns-in-readlengthdistribution#latest


## Input

input is a set of bam file or a file with suffix '.list' containing the path to the bam


## Example

```
$ java -jar dist/samreadlengthdistribution.jar -m INSERT_LENGTH src/test/resources/S*.bam 
#sample	COUNT-INSERT_LENGTH	MIN-INSERT_LENGTH	MAX-INSERT_LENGTH	AVG-INSERT_LENGTH	MEDIAN-INSERT_LENGTH	[-Inf/0[	[0/10[	[10/20[	[20/30[	[30/40[	[40/50[	[50/100[	[100/150[	[150/200[	[200/250[[250/300[	[300/350[	[350/400[	[400/450[	[450/500[	[500/1000[	[1000/Inf[
S1	999	310	696	501.31	502.00	0	0	0	0	0	0	0	0	0	003	18	115	336	527	0
S2	999	322	667	500.69	501.00	0	0	0	0	0	0	0	0	0	002	17	138	332	510	0
S3	999	322	667	500.69	501.00	0	0	0	0	0	0	0	0	0	002	17	138	332	510	0
S4	999	368	654	500.76	501.00	0	0	0	0	0	0	0	0	0	000	23	138	327	511	0
S5	999	333	641	499.78	501.00	0	0	0	0	0	0	0	0	0	004	26	126	330	513	0
```


