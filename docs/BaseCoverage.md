# BaseCoverage

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

'Depth of Coverage' per base.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar basecoverage  [options] Files

Usage: basecoverage [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -Z, --disable-normalization
      disable depth normalization _on median
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, --region, --interval
      An interval as the following syntax : "chrom:start-end". Some jvarkit 
      programs also allow the following syntax : "chrom:middle+extend"  or 
      "chrom:start-end+extend" or "chrom:start-end+extend-percent%".A program 
      might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)
    --mapq
       min mapping quality.
      Default: 1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --runmed
       moving median size
      Default: 31
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * depth
 * bam
 * sam
 * coverage
 * vcf



## Creation Date

20220420

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/basecoverage/BaseCoverage.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/basecoverage/BaseCoverage.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **basecoverage** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of bam files or one file with the '.list' suffix containing the path to the bams

## Example:

```
$ java -jar dist/basecoverage.jar --region "RF03:491-500" -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 

lindenb@asimov:~/src/jvarkit$ java -jar dist/basecoverage.jar --bed jeter.bed -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 

##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##basecoverage.meta=compilation:20220420110509 githash:afbe74ab2 htsjdk:2.24.1 date:20220420110554 cmd:--bed jeter.bed -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF03	491	.	N	.	.	.	DP=28	DP	3	6	6	8	5
RF03	492	.	N	.	.	.	DP=27	DP	3	5	5	8	6
RF03	493	.	N	.	.	.	DP=27	DP	3	5	5	8	6
RF03	494	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	495	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	496	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	497	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	498	.	N	.	.	.	DP=24	DP	3	4	4	8	5
RF03	499	.	N	.	.	.	DP=25	DP	3	4	4	9	5
RF03	500	.	N	.	.	.	DP=23	DP	3	3	3	9	5


```


