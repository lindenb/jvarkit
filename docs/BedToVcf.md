# BedToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bed2vcf  [options] Files

Usage: bed2vcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    -s, --symbolic
      use symbolic REF allele if bed length > 'x'. 'x' <=0 to disable
      Default: 100
    --version
      print version and exit

```


## Keywords

 * bed
 * vcf



## Creation Date

20240604

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2vcf/BedToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2vcf/BedToVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bed2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele.
The VCF can be used later for genotyping.

input is a set of BED file or one file with a '.list' extension containing the path to the BED files.

## Example



```
$ bcftools query -f '%CHROM\t%POS0\t%END\n' src/test/resources/rotavirus_rf.vcf.gz |\
	java -jar dist/jvarkit.jar bed2vcf -R src/test/resources/rotavirus_rf.fa | head -n 30

##fileformat=VCFv4.2
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
RF01	970	.	A	N	.	.	.
RF02	251	.	A	N	.	.	.
RF02	578	.	G	N	.	.	.
RF02	877	.	T	N	.	.	.
RF02	1726	.	T	N	.	.	.
RF02	1962	.	TACA	N	.	.	.
RF03	1221	.	C	N	.	.	.
RF03	1242	.	C	N	.	.	.
RF03	1688	.	T	N	.	.	.
RF03	1708	.	G	N	.	.	.
RF03	2150	.	T	N	.	.	.
RF03	2201	.	G	N	.	.	.
RF03	2315	.	G	N	.	.	.
RF03	2573	.	A	N	.	.	.
RF04	887	.	A	N	.	.	.
RF04	991	.	T	N	.	.	.
RF04	1241	.	T	N	.	.	.
```


