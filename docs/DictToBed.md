# DictToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert a SAM dictionary from vcf,sam,bam,dict, etc.. to bed.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar dict2bed  [options] Files

Usage: dict2bed [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-errors
      ignore errors, skip files that don't have a dictionary
      Default: false
    --no-header
      disable print header
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --skip-attributes
      skip attributes output
      Default: false
    --version
      print version and exit

```


## Keywords

 * dict
 * bed
 * sam
 * bam
 * vcf



## Creation Date

20240603

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2bed/DictToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2bed/DictToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **dict2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to bed

## Example

```
$ java -jar dist/jvarkit.jar dict2bed ~/src/jvarkit-git/src/test/resources/*.bam |  head | column -t
chrom  start  end        path                                                                         tid  buildName  AS  M5  SP  UR
chr1   0      248956422  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  0    GRCh38     .   .   .   .
chr2   0      242193529  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  1    GRCh38     .   .   .   .
chr3   0      198295559  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  2    GRCh38     .   .   .   .
chr4   0      190214555  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  3    GRCh38     .   .   .   .
chr5   0      181538259  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  4    GRCh38     .   .   .   .
chr6   0      170805979  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  5    GRCh38     .   .   .   .
chr7   0      159345973  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  6    GRCh38     .   .   .   .
chr8   0      145138636  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  7    GRCh38     .   .   .   .
chr9   0      138394717  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  8    GRCh38     .   .   .   .
```

```
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2bed | head | column -t
chrom  start  end        path                                                     tid  buildName
1      0      249250621  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  0    GRCh37
2      0      243199373  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  1    GRCh37
3      0      198022430  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  2    GRCh37
4      0      191154276  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  3    GRCh37
5      0      180915260  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  4    GRCh37
6      0      171115067  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  5    GRCh37
7      0      159138663  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  6    GRCh37
8      0      146364022  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  7    GRCh37
9      0      141213431  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  8    GRCh37
```


