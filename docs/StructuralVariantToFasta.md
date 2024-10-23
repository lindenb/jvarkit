# StructuralVariantToFasta

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert VCF of structural variant(s) to fasta for pggb


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar sv2fasta  [options] Files

Usage: sv2fasta [options] Files
  Options:
    --exclude-filtered
      Exclude FILTER-ed variants.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, --region, --interval
      interval. An interval as the following syntax : "chrom:start-end". Some 
      jvarkit programs also allow the following syntax : "chrom:middle+extend" 
       or "chrom:start-end+extend" or "chrom:start-end+extend-percent%".A 
      program might use a Reference sequence to fix the chromosome name (e.g: 
      1->chr1) 
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --skip-no-sv
      Exclude VCF without structural variant
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * cnv
 * fasta



## Creation Date

20230403

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sv2fasta/StructuralVariantToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sv2fasta/StructuralVariantToFasta.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sv2fasta** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# motivation

create a fasta input for `pggb` (pangenome graph builder)

# input

input is a set of indexed VCF files or a file with the suffix '.list' containing the path to the VCFs

# example

```
java -jar dist/jvarkit.jar sv2fasta --interval "chr1:204091-325320" --skip-no-sv --exclude-filtered -R ref.fasta vcfs.list > tmp.fa

samtools faidx tmp.fa

pggb -i tmp.fa \
	-o OUT \
    -n `grep ">" tmp.fa | wc -l` \
	-p 90 -s 100 \
	-t 5 1>&2

```



