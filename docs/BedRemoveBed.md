# BedRemoveBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Remove bed file from each record of input bed file. Output is a SETFILE


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedremovebed  [options] Files

Usage: bedremovebed [options] Files
  Options:
    -exclude, -X, --exclude-bed
      Bed to be substracted from original input (gene)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -input, -I, --include-bed
      Bed to be keep from original input (gene)
    -o, --out
      Output is a setfile. Output file. Optional . Default: stdout
    -R, -r, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --version
      print version and exit

```


## Keywords

 * bed



## Creation Date

20221210

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedtools/BedRemoveBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedtools/BedRemoveBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedremovebed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

remove bed from each record of bed

## EXAMPLES

### Example



