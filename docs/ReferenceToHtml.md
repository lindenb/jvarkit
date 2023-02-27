# ReferenceToHtml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Produce a html view of eah interval from a VCF file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar ref2html  [options] Files

Usage: ref2html [options] Files
  Options:
    ---prefix
      file prefix.
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, --region, --regions, --interval
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
  * -o, --output
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * vcf
 * html



## Creation Date

20211012

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ref2html/ReferenceToHtml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ref2html/ReferenceToHtml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ref2html** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


