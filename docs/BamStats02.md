# BamStats02

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Statistics about the flags and reads in a BAM


## Usage

```
Usage: bamstats02 [options] Files
  Options:
    -B, --bed
      Optional Bed File
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamstats02
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats02.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats02.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamstats02** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$  find dir -name "*final.bam" | xargs  java -jar dist/bamstats02.jar -B capture.bed  > output.tsv
$  verticalize output.tsv

>>> 2
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      2
$4      mapq    60
$5      inTarget        1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      1
$8      READ_UNMAPPED   0
$9      MATE_UNMAPPED   0
$10     READ_REVERSE_STRAND     1
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   0
$13     SECOND_IN_PAIR  1
$14     NOT_PRIMARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   463982
<<< 2

>>> 3


>>> 3
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      .
$4      mapq    0
$5      inTarget        -1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      0
$8      READ_UNMAPPED   1
$9      MATE_UNMAPPED   1
$10     READ_REVERSE_STRAND     0
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   1
$13     SECOND_IN_PAIR  0
$14     NOT_PRIMARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   458630
<<< 3

>>> 4
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      .
$4      mapq    0
$5      inTarget        -1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      0
$8      READ_UNMAPPED   1
$9      MATE_UNMAPPED   1
$10     READ_REVERSE_STRAND     0
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   0
$13     SECOND_IN_PAIR  1
$14     SECONDARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   458630
<<< 4
```
```





### See also

BamStats02View


