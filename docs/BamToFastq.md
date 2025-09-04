# BamToFastq

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert paired-end SAM to fastq using a memory buffer.


## DEPRECATED

Use samtools collate | samtools fastq

## Usage

```
Usage: java -jar dist/bam2fastq.jar  [options] Files
Usage: bam2fastq [options] Files
  Options:
    -d, --distance
      put the reads in memory if they're lying within that distance. A 
      distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 5000
    -R1, --forward
      Save fastq_R1 to file (default: stdout)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    -R2, --reverse
      Save fastq_R2 to file (default: interlaced with forward)
    -R0, --single
      Save single-end to this file. If unspecified, single-end reads are 
      ignored. 
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -U1, --unpaired-forward
      Save unresolved forward pair to file. If unspecified, unresolved reads 
      are ignored.
    -U2, --unpaired-reverse
      Save unresolved forward pair to file. If unspecified, unresolved reads 
      are ignored.
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * fastq
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2fastq
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20131120

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2fastq** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## DEPRECATED

use  
```
 samtools collate | samtools fastq
```

## Example

```
$ java -jar dist/bam2fastq.jar  src/test/resources/S2.bam  | head -n 16
@RF01_38_466_2:0:0_3:0:0_8f
TCTAGTCAGAATATTTATCATTTATATATAACTCACAATCCGCATTTCAAATTCCAATATACTATTCTTC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@RF01_38_466_2:0:0_3:0:0_8f
CCAACCAGAACATAACTGCATTTAAATTTGATGATAATTAAGTTAAACTTGCTGGATCCATCAATTAATC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@RF01_20_472_1:0:0_3:0:0_32
GTTTTACCCACCAGAACATAACTGCATTTAAATTTGATGATAATGAAGTTAAAATTGCTGGCTCCATCAA
+
2222222222222222222222222222222222222222222222222222222222222222222222
@RF01_20_472_1:0:0_3:0:0_32
TGGGGAAGTATAATCTAATCTTGTCAGAATATTTATCATTTATATATAACTCACAATCCGCAGTTCAACT
+
2222222222222222222222222222222222222222222222222222222222222222222222
```


## Cited In:

  * "Plastomes of nine hornbeams and phylogenetic implications", Ying Li & al;  Ecology and Evolution, 2018; DOI: 10.1002/ece3.4414; https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.4414 


