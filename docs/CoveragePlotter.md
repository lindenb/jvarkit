# CoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find anomaly of depth in intervals+bams


## Usage

```
Usage: coverageplotter [options] Files
  Options:
    --alpha
      line opacity. A decimal number between 0.0 and 1.0. If the value ends 
      with '%' it is interpretted as a percentage eg. '1%' => '0.01'. A slash 
      '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 1.0
    --arc-alpha
      arc opacity. A decimal number between 0.0 and 1.0. If the value ends 
      with '%' it is interpretted as a percentage eg. '1%' => '0.01'. A slash 
      '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 1.0
  * -B, --bams
      list of bams. one file with a '.list' suffix is interpretted a a list of 
      path to the bams
      Default: []
    --dimension
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or it can be the path to an existing png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=300]
    --black, --exclude
      Optional. BED Tabix indexed black-listed region
    --extend, -x
      extend original interval by this fraction
      Default: 1.0
    --gtf
      Optional Tabix indexed GTF file. Will be used to retrieve an interval by 
      gene name, or to display gene names in a region.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-known-containing
      Ignore known CNV containing the whole region (prevent large known CNV to 
      be displayed)
      Default: false
    --known
      Optional Tabix indexed Bed or VCF file containing known CNV. Both types 
      must be indexed.
    --manifest
      Optional. Manifest file
    --mapq
      min mapping quality
      Default: 1
    --max-arc
      max arc length in bp.
      Default: 10000000
    --max-depth
      ignore position if depth > 'x'
      Default: 500
    --min-arc
      min arc length in bp.
      Default: 1000
    -o, --output
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    --prefix
      Image File Prefix.
      Default: <empty string>
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --rrff
      Only display arcs where the strands the the read and its mate are 
      Forward-Forward or Reverse-Reverse
      Default: false
    --skip-center
      When calculating the median depth, only consider the extended region, 
      not the original interval.
      Default: false
    --smooth
      sliding window smooth size.
      Default: 250
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * depth
 * coverage


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

input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)


```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```


