# BamToSVG

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

BAM to Scalar Vector Graphics (SVG)


## Usage

```
Usage: bam2svg [options] Files
  Options:
    --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: Accept All/ Filter out nothing
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -i, --interval, --region
      interval CHROM:START-END
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -c, --showclipping
      Show clipping
      Default: false
    -S, --vcf
      add VCF indexed with tabix. Optinal. the Samples's name must be the same 
      than in the BAM
      Default: []
    --version
      print version and exit
    -w, --width
      Page width
      Default: 1000

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * svg


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2svg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/BamToSVG.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/BamToSVG.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/BamToSVGTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/BamToSVGTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ java -jar dist/bam2svg.jar \
    -R human_g1k_v37.fasta \
    -i "19:252-260" \
    -S variants.vcf.gz \
    file.bam > out.svg
```

## Gallery

https://twitter.com/yokofakun/status/523031098541232128

![bam2svg](https://pbs.twimg.com/media/B0IuAw2IgAAYfNM.jpg)

https://twitter.com/yokofakun/status/522415314425090048

![bam2svg-2](https://pbs.twimg.com/media/Bz_99ayIMAAK57s.jpg)



