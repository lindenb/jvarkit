# BamStats01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Statistics about the reads in a BAM.


## Usage

```
Usage: samstats01 [options] Files
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
    -o, --output
      Output file. Optional . Default: stdout
    -q, --qual
      min mapping quality
      Default: 30
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary 
    -B, --bed, --regions
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
    --version
      print version and exit

```


## Keywords

 * sam
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samstats01
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20130705

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samstats01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$ java -jar dist/samstats01.jar --bed "RF03:1-3000"  src/test/resources/S*.bam 
#Filename	Sample	UNMAPPED	NOT_PRIMARY	FAIL_VENDOR_QUALITY	OFF_TARGET	DUPLICATE	FAIL_MAPPING_QUALITY	OK_CALLING
src/test/resources/S1.bam	S1	0	0	0	1718	0	0280
src/test/resources/S2.bam	S2	0	0	0	1718	0	0280
src/test/resources/S3.bam	S3	0	0	0	1718	0	0280
src/test/resources/S4.bam	S4	0	0	0	1718	0	0280
src/test/resources/S5.bam	S5	0	0	0	1718	0	0280
```


