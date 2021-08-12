# WGSCoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Whole genome coverage plotter


## Usage

```
Usage: wgscoverageplotter [options] Files
  Options:
    --clip, --cap
      Don't show coverage to be greater than 'max-depth' in the SVG file.
      Default: false
    --dimension
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or it can be the path to an existing png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=500]
    --disable-paired-overlap
      Count overlapping bases with mate for paired-end
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -I, --include-contig-regex
      Only keep chromosomes matching this regular expression. Ignore if blank.
      Default: <empty string>
    --mapq
      min mapping quality
      Default: 1
    -C, --max-depth
      Max depth to display. The special value '-1' will first compute the 
      average depth and the set the max depth to 2*average
      Default: 100
    --min-contig-length
      Skip chromosome with length < 'x'. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      When using the option --samples, use this partition Data partitioning 
      using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    --percentile
      How to we bin the coverage under one pixel.
      Default: median
      Possible Values: [median, average, min, max]
    --points
      Plot the coverage using points instead of areas.
      Default: false
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --samples
      Limit to those groups. See also --partition. Multiple separated with 
      commas. 
      Default: <empty string>
    -X, --skip-contig-regex
      Skip chromosomes matching this regular expression. Ignore if blank.
      Default: <empty string>
    --version
      print version and exit
    -D
      set some css style elements. '-Dkey=value'. Undocumented.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * svg
 * bam
 * depth
 * coverage



## See also in Biostars

 * [https://www.biostars.org/p/104063](https://www.biostars.org/p/104063)
 * [https://www.biostars.org/p/475162](https://www.biostars.org/p/475162)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew wgscoverageplotter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20201125

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **wgscoverageplotter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Files

Input is an indexed BAM or CRAM file

Output is a SVG file

## Example
```
java -jar dist/wgscoverageplotter.jar --dimension 1500x500 -C -1 --clip -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam --include-contig-regex "RF.*" --percentile median  > ~/jeter.svg
```

## Screenshot

https://twitter.com/yokofakun/status/1331898068002861056

![twitter](https://pbs.twimg.com/media/EnvaOnNW4AAkGTz?format=jpg&name=medium "Screenshot")

