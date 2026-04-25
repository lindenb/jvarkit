# BedToHilbert

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

BED to hilbert curve as SVG.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bed2hilbert  [options] Files

Usage: bed2hilbert [options] Files
  Options:
    --columns
      Optional Comma separated list names of column starting from the 3rd 
      column (after 'end'). Use '.' to ignore.
      Default: <empty string>
    --css
      load custom CSS style file
    --first-line-header
      The first line of the bed file is a header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --min-contig-length
      keep chromosomes which length is greater than 'x'
      Default: 0
    --omit-xml-desclaration
      Don't print XM declaration.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -r, --recursion
      Hilbert Curve level of recursion
      Default: 6
  * -R, --reference, --dict
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --regex
      keep chromosomes matching that regular expression
    --version
      print version and exit
    -w, --width
      Image width
      Default: 1000

```


## Keywords

 * bed
 * xml
 * svg
 * hilbert



## Creation Date

20260417

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2hilbert/BedToHilbert.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2hilbert/BedToHilbert.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2hilbert/BedToHilbertTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2hilbert/BedToHilbertTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bed2hilbert** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Creates a Hilbert SVG Graph  for a BED file.


## Example

```
 java -jar dist/jvarkit.jar bed2hilbert  -R src/test/resources/rotavirus_rf.dict input.bed > out.svg
```


