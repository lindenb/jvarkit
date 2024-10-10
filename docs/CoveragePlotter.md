# CoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Display an image of depth to display any anomaly an intervals+bams


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar coverageplotter  [options] Files

Usage: coverageplotter [options] Files
  Options:
    --css
      Custom CSS file. format <sample> <css>. One per line. eg. "sample1 
      stroke:red; 
    --dimension, --dim
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or the word 'screen' or it can be the path to an existing 
      png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=300]
    --extend, -x
      Extending interval. The following syntaxes are supported: 1000; 1kb; 
      1,000; 30%(shrink); 150% (extend); 0.5 (shrink); 1.5 (extend)
      Default: 3.0
    --gff, --gff3
      Optional Tabix indexed GFF3 file. Will be used to retrieve an interval 
      by gene name, or to display gene names in a region.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-known-containing
      Ignore known CNV containing the whole region (prevent large known CNV to 
      be displayed)
      Default: false
    --include-center
      When calculating the median depth, also consider the original user's 
      region, not only the extended interval.
      Default: false
  * --regions, --region, --interval
      Interval region
    --known
      Optional Tabix indexed Bed or VCF file containing known CNV. Both types 
      must be indexed.
    --loess
      Run Loess smoothing on GC%. Experimental. For now, I find the smooting 
      is too strong.
      Default: false
    --mapq
      min mapping quality
      Default: 1
    --max-depth
      ignore position if depth > 'x'
      Default: 500
    --max-y
      Max normalized Y
      Default: 3.0
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --smooth
      Run median smooth on this number of pixels. (ignore if <=1)
      Default: 10
    --svg-only
      Force SVG-only output (default is HTML+SVG).
      Default: false
    --use-average
      Calculating the median depth can be memory consumming for large regions. 
      If the region is larger than 'x', use 'average' instead of 'median'. A 
      distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 2000000
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * depth
 * coverage
 * svg



## See also in Biostars

 * [https://www.biostars.org/p/9536274](https://www.biostars.org/p/9536274)



## Creation Date

20200605

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/CoveragePlotter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/CoveragePlotter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/CoveragePlotterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/CoveragePlotterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **coverageplotter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

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
java -jar dist/jvarkit.jar coverageplotter -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" in.list
```

## Screenshot

!(https://pbs.twimg.com/media/Fac3XR3aAAEJoXu?format=jpg&name=medium)[https://twitter.com/yokofakun/status/1560276675614887937]



