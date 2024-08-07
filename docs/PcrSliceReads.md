# PcrSliceReads

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Mark PCR reads to their PCR amplicon


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar pcrslicereads  [options] Files

Usage: pcrslicereads [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -B, --bed
      bed file containing non-overlapping PCR fragments. Column name is 
      required. 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --random
       random seed (-1 == timestamp)
      Default: -1
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit
    -a
       if a read is mapped on multiple PCR fragments, how to resolve ambiguity
      Default: closer
      Possible Values: [zero, random, closer]
    -c
      clip read to their PCR fragments. see 
      https://github.com/lindenb/jvarkit/wiki/PcrClipReads 
      Default: false

```


## Keywords

 * pcr
 * sam
 * bam
 * cigar



## See also in Biostars

 * [https://www.biostars.org/p/149687](https://www.biostars.org/p/149687)



## Creation Date

20150707

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrSliceReads.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrSliceReads.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pcrslicereads** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## See also

* (2021) `samtools ampliconclip` - clip reads using a BED file  http://www.htslib.org/doc/samtools-ampliconclip.html

## Motivation

Mark PCR reads to their PCR amplicon https://www.biostars.org/p/149687/

Experimental, use with care.

* reads must have a read group id (see picard  addorreplacereadgroups http://broadinstitute.github.io/picard/ )
* if the read is not a paired-end read , then the quality of the read is set to zero
* if the MATE is not mapped , then the quality of the current read is set to zero
* if the properly-paired flag is not set, then the quality of the current read is set to zero
* if the read and the mate are not mapped on the same chromosome, then the quality of the current read is set to zero
* if the read and the mate are mapped on the same strand, then the quality of the current read is set to zero
* if the read(POS) < mate(POS) and read is mapped on the negative strand, then the quality of the current read is set to zero
* if the read(POS) > mate(POS) and read is mapped on the positive strand, then the quality of the current read is set to zero
* if no BED fragment is found overlapping the read, then the quality of the read is set to zero
* if a PCR fragment is found, the read group 'ID' is changed to 'ID_fragmentname'


## Example

create a bed file with bedtools, add a column 'PCR%%' as the name of the PCR fragment. Add a read group to ex1.bam , change the readgroup with pcrslicereads, select mapped reads having a MAPQ>0, sort and index, get the coverage with GATK/DepthOfCoverage using **--partitionType readgroup **


```makefile
.SHELL=/bin/bash


coverage : clipped.bam
	java -jar gatk/3.3.0/GenomeAnalysisTK.jar \
	   -T DepthOfCoverage   -omitBaseOutput \
	   -R samtools-0.1.19/examples/ex1.fa \
	   --partitionType readgroup -I $< -o coverage
	head  coverage.read_group_summary
	tail coverage.read_group_summary

clipped.bam : dist/pcrslicereads.jar test.bed
	java -jar picard/picard-tools-1.126/picard.jar  AddOrReplaceReadGroups \
		I=samtools-0.1.19/examples/ex1.bam \
		O=test.bam ID=myid LB=mylb PL=illumina PU=mypu SM=mysample
	java -jar $< -c   -B  test.bed -b test.bam |\
	samtools  view -q 1 -F 4 -bu  -  |\
	samtools  sort - clipped && samtools index clipped.bam

dist/pcrslicereads.jar:   src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrSliceReads.java
	make pcrslicereads

test.bed: samtools-0.1.19/examples/ex1.fa.fai
	bedtools makewindows -g $< -w 200 -s 50 |\
		awk '{printf("%s\tPCR%d\n",$$0,NR);}' > $@


```

```
$ head test.bed 
seq1    0       200     PCR1
seq1    50      250     PCR2
seq1    100     300     PCR3
seq1    150     350     PCR4
seq1    200     400     PCR5
seq1    250     450     PCR6
seq1    300     500     PCR7
seq1    350     550     PCR8
seq1    400     600     PCR9
seq1    450     650     PCR10
```

```
$ samtools view -h clipped.bam | more
@HD     VN:1.4  SO:coordinate
@SQ     SN:seq1 LN:1575 UR:file:/commun/data/packages/samtools-0.1.19/examples/ex1.fa   AS:ex1  M5:426e31835a6dfdcbf6c534671edf02f7
@SQ     SN:seq2 LN:1584 UR:file:/commun/data/packages/samtools-0.1.19/examples/ex1.fa   AS:ex1  M5:b6853ffe730ece50076db834dea18e3b
@RG     ID:myid PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR1    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR2    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR3    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR4    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR5    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR6    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR7    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR8    PL:illumina     PU:mypu LB:mylb SM:mysample
@RG     ID:myid_PCR9    PL:illumina     PU:mypu LB:mylb SM:mysample
(...)
EAS54_71:8:105:854:975  83      seq2    1523    71      28M5S   =       1354    -202    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTG       <<<<;<:<<;<&<;<<<<<<<<<<<<<<<<<<<       H0:i:85 H1:i:85 MF:i:18 R
G:Z:myid_PCR60  NM:i:0  UQ:i:0  Aq:i:0
EAS139_11:7:50:1229:1313        83      seq2    1528    77      23M12S  =       1376    -187    TTTTTTCTTTTTTTTTTTTTTTTTTTTGCATGCCA     <<<<,<&<7<<<<<<<<<<<<<<<<<<<<<<<<<<     H0:i:3  H1:i:7  M
F:i:18  RG:Z:myid_PCR60 NM:i:1  UQ:i:11 Aq:i:0
EAS54_65:3:320:20:250   147     seq2    1532    77      19M16S  =       1367    -200    TTTTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAA     +'''/<<<<7:;+<;::<<<;;<<<<<<<<<<<<<     H0:i:1  H1:i:2  MF:i:18 R
G:Z:myid_PCR60  NM:i:2  UQ:i:24 Aq:i:6
EAS114_26:7:37:79:581   83      seq2    1533    68      18M17S  =       1349    -219    TTTTTTTTTTTTTTTTTTTTTTTCATGCCAGAAAA     3,,,===6===<===<;=====-============     H0:i:0  H1:i:1  MF:i:18 R
G:Z:myid_PCR60  NM:i:2  UQ:i:23 Aq:i:27
```

```
$ head coverage.read_group_summary
sample_id       total   mean    granular_third_quartile granular_median granular_first_quartile %_bases_above_15
mysample_rg_myid_PCR2   1106    0.35    1       1       1       0.2
mysample_rg_myid_PCR3   674     0.21    1       1       1       0.0
mysample_rg_myid_PCR1   51      0.02    1       1       1       0.0
mysample_rg_myid_PCR6   1279    0.40    1       1       1       0.1
mysample_rg_myid_PCR7   1554    0.49    1       1       1       1.3
mysample_rg_myid_PCR4   1147    0.36    1       1       1       0.8
mysample_rg_myid_PCR5   1738    0.55    1       1       1       1.5
mysample_rg_myid_PCR9   1260    0.40    1       1       1       0.8
mysample_rg_myid_PCR8   1930    0.61    1       1       1       2.1
```



