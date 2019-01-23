# Bam2Wig

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Bam to fixedStep Wiggle converter , or BED GRAPH. Parses the cigar String to get the depth. Memory intensive: must alloc sizeof(int)*size(chrom)


## Usage

```
Usage: bam2wig [options] Files
  Options:
    -bg, --bedgraph
      Produce a BED GRAPH instead of a WIGGLE file.
      Default: false
    --display
      What kind of data should we display ?
      Default: COVERAGE
      Possible Values: [COVERAGE, CLIPPING, INSERTION, DELETION, READ_GROUPS, CASE_CTRL]
    --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: Accept All/ Filter out nothing
    -f, --format
      `Printf` Format for the values. see 
      https://docs.oracle.com/javase/tutorial/java/data/numberformat.html . 
      Use "%.01f" to print an integer. "%e" for scientific notation.
      Default: %.3f
    -t, --header
      print a UCSC custom track header: something lile track type=track_type 
      name="__REPLACE_WIG_NAME__" description="__REPLACE_WIG_DESC__". Use 
      `sed` to replace the tokens. e.g: `sed 
      '/^track/s/__REPLACE_WIG_NAME__/My data/'`
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --region, --interval
      Limit analysis to this interval. An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or "chrom:start-end+extend" 
      or "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
    --mindepth, --mindp
      When using display READ_GROUPS, What is the minimal read depth that 
      should be considered ?
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      When using display READ_GROUPS, how should we partition the ReadGroup ? 
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    --pedigree, -ped
      Pedigree file for CASE_CTRL. A pedigree is a text file delimited with 
      tabs. No header. Columns are (1) Family (2) Individual-ID (3) Father Id 
      or '0' (4) Mother Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) 
      Status : 0 unaffected, 1 affected,-9 unknown
    --percentile
      How to group data in the sliding window ?
      Default: AVERAGE
      Possible Values: [MIN, MAX, MEDIAN, AVERAGE, RANDOM, SUM]
    --version
      print version and exit
    -s, --windowShift
      window shift
      Default: 25
    -w, --windowSize
      window size
      Default: 100

```


## Keywords

 * bam
 * wig
 * wiggle
 * bed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2wig
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2wig/Bam2Wig.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2wig/Bam2Wig.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2wig** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is stdin; or  one or more BAM file sorted on coordinate; or a file ending with '.list' and containing the PATH to some bam files.

About Wiggle: [https://genome.ucsc.edu/goldenpath/help/wiggle.html](https://genome.ucsc.edu/goldenpath/help/wiggle.html)

About BedGraph: [https://genome.ucsc.edu/goldenpath/help/bedgraph.html](https://genome.ucsc.edu/goldenpath/help/bedgraph.html)

## Memory

warning: the program is memory consuming, it allocates on array of integer of the size of your longest contig.

## History:

20171115: removed cast_to_integer replaced by 'format', added percentile. Removed options --zerolength and --mindepth.

## Aggregators:

* COVERAGE :  coverage, all sample merged
* CLIPPING : consider only clipped base
* INSERTION: consider only Cigar events I. Only one base in the reference is flagged
* DELETION: consider Cigar events 'N' and 'D':
* READ_GROUPS : Number of 'samples' having a depth greater than 'min-depth'
* CASE_CTRL: ratio median(coverage-cases)/median(coverage-controls)

## Example
the input file

```bash
java -jar dist/bam2wig.jar -w 1 -s 3   examples/toy.bam
```

```
track type=wiggle_0 name="__REPLACE_WIG_NAME__" description="__REPLACE_WIG_DESC__"
fixedStep chrom=ref start=7 step=3 span=1
1
3
3
3
1
1
0
0
1
0
2
2
1
fixedStep chrom=ref2 start=1 step=3 span=1
1
2
3
4
5
6
6
5
4
3
3
```

