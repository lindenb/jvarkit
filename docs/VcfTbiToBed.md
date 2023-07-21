# VcfTbiToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

extracts BED for each contig in a tabix-indexed VCF peeking first of last variant for each chromosome.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcftbi2bed  [options] Files

Usage: vcftbi2bed [options] Files
  Options:
    --POS, -POS
      for the end position, default is to use the END position (think about 
      SV, INDEL...) of the variant. This option just use the POS.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * vcf
 * tabix



## Creation Date

20230214

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tbi2bed/VcfTbiToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tbi2bed/VcfTbiToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcftbi2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Extract the first and last variants of a tabix-index-VCF file for each chromosome.
Output is a BED file contig/start/end/vcf.
Input can be:

  - stdin (the path to the vcfs, one per line)
  - the vcfs
  - a file with the suffix '.list' containing the path to the vcfs, one per line.

## Example:



