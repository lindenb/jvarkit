# CoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Display an image of depth to display any anomaly an intervals+bams


## Usage

```
Usage: java -jar dist/coverageplotter.jar  [options] Files
Usage: coverageplotter [options] Files
  Options:
    --css
      Custom CSS file. format <sample> <css>. One per line. eg. "sample1 
      {stroke:red;} 
    --dimension, --dim
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or it can be the path to an existing png,jpg,xcf,svg file.
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
  * --region, --interval
      Interval region
    --known
      Optional Tabix indexed Bed or VCF file containing known CNV. Both types 
      must be indexed.
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
      faidx and with picard CreateSequenceDictionary
    --svg-only
      Force SVG-only output (default is HTML+SVG).
      Default: false
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * depth
 * coverage
 * svg


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew coverageplotter
```

The java jar file will be installed in the `dist` directory.


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

input a set of bam/cram files or one file with the suffix '.list' containing the path to the bams

output is a HTML+SVG file

```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" src/test/resources/*.bam 
```


