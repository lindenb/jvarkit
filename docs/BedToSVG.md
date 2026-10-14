# BedToSVG

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert BED file(s) to SVG


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bed2svg  [options] Files

Usage: bed2svg [options] Files
  Options:
    --contig-min-length
      just show chromosomes having a length greater than 'x'
      Default: 0
    --contig-regex
      just show chromosomes matching this regex
      Default: (chr)?[0-9XY]+
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --jexl-filter
      ]A JEXL Expression returning a string JEXL stands for Java EXpression 
      Language.  See 
      https://commons.apache.org/proper/commons-jexl/reference/syntax.html . 
      The variable 'record' is the current observed BedLine (or Variant) ).  
      The variable 'filename' is the current filename
      Default: <empty string>
    --max-feature-length
      discard  features where the length of reference is lower than x. A 
      distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    --min-feature-length
      discard  features where the length of reference is greater than x. A 
      distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 0
    -o, --out
      Output file. Optional . Default: stdout
    --regions
      regions of interest, multiple comma separated.
      Default: []
    --regions-file
      regions of interest, multiple comma separated.
    --version
      print version and exit
    --width
      window width in pixels
      Default: 1000
  * -R
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)

```


## Keywords

 * bed
 * svg



## Creation Date

20261013

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2svg/BedToSVG.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2svg/BedToSVG.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2svg/BedToSVGTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2svg/BedToSVGTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bed2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Convert BED file to SVG

## Example

```

```


