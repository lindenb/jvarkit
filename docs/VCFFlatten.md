# VCFFlatten

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Flatten variants to one variant


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfflatten  [options] Files

Usage: vcfflatten [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --gene-extractor
      Activate default gene extractors. Variant will be grouped by gene using 
      snpeff/bcftools/vep annotations
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --id
      Default Variant ID
      Default: FLATTEN_VARIANT
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * contrast



## Creation Date

20230222

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfflatten/VCFFlatten.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfflatten/VCFFlatten.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfflatten** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

flatten all variants to one variant. At the end, a Sample will be set to 1/1 if one or more of his Genotypes
in the VCF is contains an NON-(REF|NO_CALL) allele
The idea is to use the VCF output of this tool
in order to pipe it into `bcftools contrast`: we'll get a p-value for the whole VCF


# Example

```
$ java -jar dist/jvarkit.jar vcfflatten --gene-extractor src/test/resources/rotavirus_rf.ann.vcf.gz 2> /dev/null 
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=MULTIPLE_CONTIG,Number=0,Type=Flag,Description="Record spans multiple chromosomes. Only first chromosome is reported">
##INFO=<ID=N_VARIANTS,Number=1,Type=Integer,Description="Number of variants">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	970	.	N	<FLATTEN_VARIANT>	.	.	MULTIPLE_CONTIG;N_VARIANTS=45	GT	1/1	1/1	1/1	1/1	1/1
RF01	970	.	N	<ANN/GeneName:Gene_18_3284>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF01	970	.	N	<ANN/FeatureId:AAA47319.1>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF01	970	.	N	<ANN/GeneId:Gene_18_3284>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF02	1962	.	N	<ANN/GeneId:UniProtKB/Swiss-Prot:P12472>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/FeatureId:CAA32215.1>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/FeatureId:CAA32213.1>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/GeneName:Gene_1621_1636>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
```


