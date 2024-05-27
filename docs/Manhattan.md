# Manhattan

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Manhattan plot SVG picture from different sources.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar manhattan  [options] Files

Usage: manhattan [options] Files
  Options:
    --bed, -L
      BED file containing one or more region to observe. Use the whole 
      chromosome if undefined
    -D, --define
      Dynamic parameters
      Syntax: -Dkey=value
      Default: {}
    --dimension
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or it can be the path to an existing png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=300]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --include
      include chromosomes matching the following expression
      Default: (chr)?[0-9XY]+
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --version
      print version and exit
    -l
      List available handlers and exit
      Default: false
    -n
      handler name
      Default: default

```


## Keywords

 * chromosome
 * reference
 * chart
 * visualization
 * svg



## Creation Date

20220525

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/manhattan/Manhattan.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/manhattan/Manhattan.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **manhattan** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Node

this program is very unstable. I often change everything...



