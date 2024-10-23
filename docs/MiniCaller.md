# MiniCaller

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Simple and Stupid Variant Caller designed for @AdrienLeger2


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar minicaller  [options] Files

Usage: minicaller [options] Files
  Options:
    --bad-ad-ratio
      Filter Genotype if x< ALT/(REF+ALT) < (1-x).
      Default: 0.2
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    --gt-fraction
      ignore genotype ALT/(REF+ALT) < x
      Default: 0.05
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
      Default: 1
    --min-gt-allele-depth
      min genotype allele DP
      Default: 1
    --min-gt-depth
      min genotype DP
      Default: 1
    -d, --mindepth
      Min depth
      Default: 20
    --other-reference
      Other fasta references if you mix bam mapped on different fasta (will 
      try to convert chromosomes names). Indexed fasta Reference file. This 
      file must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Main fasta reference. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
  * -r, --region
      An interval as the following syntax : "chrom:start-end". Some jvarkit 
      programs also allow the following syntax : "chrom:middle+extend"  or 
      "chrom:start-end+extend" or "chrom:start-end+extend-percent%".A program 
      might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * calling
 * vcf



## Creation Date

201500306

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/calling/MiniCaller.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/calling/MiniCaller.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **minicaller** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Cited-In

  * "Direct Head-to-Head Evaluation of Recombinant Adeno-associated Viral Vectors Manufactured in Human versus Insect Cells". Kondratov & al. Molecular Therapy. [https://doi.org/10.1016/j.ymthe.2017.08.003](https://doi.org/10.1016/j.ymthe.2017.08.003).
  *  METHODS OF ENHANCING BIOLOGICAL POTENCY OF BACULOVIRUS SYSTEM-PRODUCED RECOMBINANT ADENO-ASSOCIATED VIRUS   United States Patent Application 20200123572 . United States Patent Application 20200123572

## Example

```bash
$  java -jar dist/minicaller.jar -R ref.fa  bam.list > out.vcf
```

```
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Depth ReforAlt|Strand : RF,RR,AF,AR">
##FORMAT=<ID=DPG,Number=G,Type=Integer,Description="Depth for each allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is indel">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	4	.	T	A	.	.	DP=65	GT:DP:DP4:DPG	0/1:22:20,0,2,0:20,2	./.	./.	./.
rotavirus	5	.	T	A	.	.	DP=83	GT:DP:DP4:DPG	0:27:27,0,0,0:27	./.	0/1:25:20,0,5,0:20,5	./.
rotavirus	6	.	T	A	.	.	DP=97	GT:DP:DP4:DPG	0:30:30,0,0,0:30	0/1:21:20,0,1,0:20,1	0/1:33:31,0,2,0:31,2	./.
rotavirus	7	.	T	A	.	.	DP=112	GT:DP:DP4:DPG	0/1:38:36,0,2,0:36,2	0/1:23:21,0,2,0:21,2	0:37:37,0,0,0:37	./.
rotavirus	8	.	A	C	.	.	DP=122	GT:DP:DP4:DPG	0/1:41:38,0,3,0:38,3	0/1:26:25,0,1,0:25,1	0/1:40:38,0,2,0:38,2	./.
rotavirus	9	.	A	C	.	.	DP=139	GT:DP:DP4:DPG	0/1:46:44,0,2,0:44,2	0/1:29:27,0,2,0:27,2	0/1:48:44,0,4,0:44,4	./.
```


