# WGSCoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Whole genome coverage plotter


## Usage

```
Usage: wgscoverageplotter [options] Files
  Options:
    --clip, --cap
      Don't allow coverage to be greater than 'max-depth' in the SVG file.
      Default: false
    --dimension
      Image Dimension.
      Default: java.awt.Dimension[width=1000,height=500]
    --disable-paired-overlap
      Disable: Count overlapping bases with mate for paired-end
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality
      Default: 1
    -C, --max-depth
      Max depth to display
      Default: 100
    --min-contig-length
      Skip chromosome with length < 'x'
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --skip-contig-regex
      Skip chromosome matching this regular expression
      Default: (NC_007605|hs37d5)
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
$ ./gradlew wgscoverageplotter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20201125

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java)


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

input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)


```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```


