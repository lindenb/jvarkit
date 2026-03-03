# MiniGenotyper

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Simple and Stupid Variant Genotyper


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar minigenotyper  [options] Files

Usage: minigenotyper [options] Files
  Options:
    --bad-ad-ratio
      Ignore Alt if x < AD ratio < (1-x)
      Default: 0.1
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --disable-overlap-detection
      Disable Paired-end read overlap detection
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -mapq, --mapq
      min mapping quality
      Default: 1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    --min-base-quality
      min base quality
      Default: -1
    --min-gt-depth
      min genotype DP
      Default: 5
    --mode
      how to scan the bam. Using random-access (ok if small number of variant) 
      or streaming (no index required, large number of variants)
      Default: random_access
      Possible Values: [random_access, streaming]
    --no-id
      don't print ID column
      Default: false
    --no-info
      don't print INFO column
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict. 
      This can be used multiple times if there is more than one REF.
      Default: []
    --skip-illegal-variant
      just skip illegal variant in the VCF (not diallelic SNPs)
      Default: false
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
  * -V, --variant
      Variants to Genotype
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * calling
 * vcf



## Creation Date

20260221

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/minigenotyper/MiniGenotyper.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/minigenotyper/MiniGenotyper.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/minigenotyper/MiniGenotyperTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/minigenotyper/MiniGenotyperTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **minigenotyper** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of indexed bam or a file with the paths ending with '.list'

## Example

```
find src/test/resources/ -type f -name "S*.bam" > jeter.list
java -jar dist/jvarkit.jar minigenotyper \
 	-V src/test/resources/rotavirus_rf.vcf.gz  \
 	-R src/test/resources/rotavirus_rf.fa \
 	jeter.list --skip-illegal-variant > out.vcf
 ```



