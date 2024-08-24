# DictToXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert a SAM dictionary from vcf,sam,bam,dict, etc.. to XML.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar dict2xml  [options] Files

Usage: dict2xml [options] Files
  Options:
    --duplicate
      keep duplicates (default behavior is to keep one dictionary
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-errors
      ignore errors, skip files that don't have a dictionary
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * dict
 * xml
 * sam
 * bam
 * vcf



## Creation Date

20240824

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2xml/DictToXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2xml/DictToXml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **dict2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to XML
Then we can use XSLT to generate code...

## Example

```
$ java -jar dist/jvarkit.jar dict2xml ~/src/jvarkit-git/src/test/resources/*.bam | xmllint --format -
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
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2xml | xmllint --format -

```

