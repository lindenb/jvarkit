# CoverageGrid

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Display an image of depth to display any anomaly an intervals+bams as a grid image


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar coveragegrid  [options] Files

Usage: coveragegrid [options] Files
  Options:
    --dimension, --dim
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or the word 'screen' or it can be the path to an existing 
      png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=300]
    --extend, -x
      extends the interval x times
      Default: 3.0
    --format
      output format
      Default: SVG
      Possible Values: [SVG, SVG_GZ, PNG, JPG, PS, PS_GZ]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * --regions, --region, --interval
      Interval region
    --mapq
      min mapping quality
      Default: 1
    --max-y
      Max normalized Y
      Default: 3.0
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --threads
      number of threads
      Default: 1
    --title
      Title
      Default: <empty string>
    --type
      Plot type
      Default: MEDIAN_COVERAGE
      Possible Values: [COVERAGE, MEDIAN_COVERAGE, PILEUP, PILEUP_PAIR]
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * depth
 * coverage
 * svg
 * postscrip



## Creation Date

20241009

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/coveragegrid/CoverageGrid.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/coveragegrid/CoverageGrid.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **coveragegrid** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## input

input a set of bam/cram files or one file with the suffix '.list' containing the path to the bams

## output

output is a HTML+SVG file

## example:

```
find dir -type f -name "*bam" > in.list 
java -jar dist/jvarkit.jar coveragegrid -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" in.list
```


