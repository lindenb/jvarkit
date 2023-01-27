# BamToMNV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

MNV haplotypes


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bam2mnv  [options] Files

Usage: bam2mnv [options] Files
  Options:
    --distance
      max distance distance between two snp.
      Default: 100
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -I, --input, --bam
      add this indexed bam
      Default: []
    --mapq
      min mapping quality
      Default: 1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --out
      Output file. Optional . Default: stdout
    --pedigree, --ped
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1|male|M:male;2|female|F:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary 
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * phased
 * genotypes
 * bam



## Creation Date

20211208

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamToMNV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamToMNV.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2mnv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ find src/test -type f -name "*.bam" > jeter.list
$ java -jar dist/bam2mnv.jar --input jeter.list src/test/resources/rotavirus_rf.vcf.gz --distance 3000

#CHROM1	POS1	REF1	ALT1	CHROM2	POS2	REF2	ALT2	distance	S1	S2	S3	S4	S5
RF03	1221	C	G	RF03	1242	C	A	22	ref	v1	v1	v2	ref
RF03	1221	C	G	RF03	1688	T	G	468	ref	v1	v1	ref	v2
RF03	1221	C	G	RF03	1708	G	T	488	ref	v1	v1	ref	v2
RF03	1242	C	A	RF03	1688	T	G	447	ref	ref	ref	v1	v2
RF03	1242	C	A	RF03	1708	G	T	467	ref	ref	ref	v1	v2
RF03	1688	T	G	RF03	1708	G	T	21	ref	ref	ref	ref	cis
RF03	1688	T	G	RF03	2150	T	A	463	ref	ref	ref	ref	v1
RF03	1688	T	G	RF03	2201	G	C	514	ref	v2	v2	.	.
RF03	1708	G	T	RF03	2150	T	A	443	.	ref	ref	ref	v1
RF03	1708	G	T	RF03	2201	G	C	494	ref	v2	v2	.	.
RF03	2150	T	A	RF03	2201	G	C	52	ref	v2	v2	ref	.
RF03	2150	T	A	RF03	2573	A	G	424	v1	v2	v2	ref	.
RF03	2201	G	C	RF03	2573	A	G	373	.	cis	cis	.	.
RF04	1900	A	C	RF04	1920	A	T	21	v2	ref	ref	ref	v1
RF05	41	T	C	RF05	499	A	T	459	ref	v2	v2	v1	ref
RF05	879	C	A	RF05	1297	T	G	419	ref	cis	cis	ref	ref
RF05	1297	T	G	RF05	1339	A	C	43	ref	cis	cis	ref	.
RF06	517	C	A	RF06	543	G	C	27	.	v2	v2	.	ref
RF06	668	A	G	RF06	695	T	C	28	ref	ref	ref	.	v1
RF07	225	C	A	RF07	684	T	G	460	ref	.	.	v2	.
RF08	926	A	C	RF08	992	G	C	67	ref	v2	v2	.	v1
RF09	294	T	C	RF09	317	C	A	24	.	ref	ref	ref	v2
```


