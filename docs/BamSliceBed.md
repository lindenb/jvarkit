# BamSliceBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

For @wouter_decoster : slice (long reads) overlapping the records of a BED file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bamslicebed  [options] Files

Usage: bamslicebed [options] Files
  Options:
    --attributes
      Leep the following attributes (separated by spaces/comma/semicolon)
      Default: <empty string>
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
  * -B, --bed, --pcr
      Regions containing non-overlapping PCR fragments. A source of intervals. 
      The following suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, 
      gff, gff.gz, gtf.gz.Otherwise it could be an empty string (no interval) 
      or a list of plain interval separated by '[ \t\n;,]'
      Default: (unspecified)
    --clip
      Do not remove the bases but soft clip them.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * bed



## Creation Date

20191030

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/BamSliceBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/BamSliceBed.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pcr/BamSliceBedTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pcr/BamSliceBedTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamslicebed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## See also

* (2021) `samtools ampliconclip` - clip reads using a BED file  http://www.htslib.org/doc/samtools-ampliconclip.html


## Input/Output

input is a bam

output is a bam


name of the BED record is appended to the original  read name/

unmapped reads, reads without cigar or reads that don't overlap any BED record are discarded

MAPQ is set to 255

reads are converted to singled end

optional args are not filled

bounding bases with cigar string without cigar operator M/X/= are discarded.



## Example

```
$ cat jeter.bed
RF01	10	15
RF01	20	25
RF01	30	35
```

```
$ java -jar dist/bamslicebed.jar -B jeter.bed ./src/test/resources/S1.bam |samtools sort -T tmp -o jeter.bam -
$ samtools view jeter.bam 
RF01_1_483_2:0:0_3:0:0_41#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:11:15	0	RF01	12	255	4M	*	0	0	CTAT	2222
RF01_1_483_2:0:0_3:0:0_41#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGC	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_1_483_2:0:0_3:0:0_41#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_27_590_3:0:0_1:0:0_68#RF01:31:35	0	RF01	31	255	5M	*	0	0	CATCT	22222

samtools tview jeter.bam src/test/resources/rotavirus_rf.fa


1         11        21        31    
ggctattaaagctatacaATGGGGAAGTATAATCTA
          .....     .....     .....
          .....     ....C     .....
          .....     .....
          .....     .....
           ....     .....
                              .....
                              .....
                              .....
                              C....
```


