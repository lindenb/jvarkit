# SamReadLengthDistribution

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sam read length distribution


## Usage

```
Usage: samreadlengthdistribution [options] Files
  Options:
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
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -w, --windows
      A 'range of integers' is a list of integers in ascending order separated 
      with semicolons.
      Default: [[-Inf/0[, [0/10[, [10/20[, [20/30[, [30/40[, [40/50[, [50/100[, [100/150[, [150/200[, [200/250[, [250/300[, [300/350[, [350/400[, [400/450[, [450/500[, [500/1000[, [1000/Inf[]
    -m
      method, how should I get the read length ?
      Default: SEQ_LENGTH
      Possible Values: [SEQ_LENGTH, CIGAR_REF_LENGTH, CIGAR_PADDED_REF_LENGTH]

```


## Keywords

 * sam
 * bam
 * histogram


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samreadlengthdistribution
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistribution.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistribution.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistributionTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamReadLengthDistributionTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samreadlengthdistribution** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Because gatk is buggy: http://gatkforums.broadinstitute.org/gatk/discussion/8342/duplicate-columns-in-readlengthdistribution#latest


## Input

input is a set of bam file or a file with suffix '.list' containing the path to the bam


## Example

```
$ java -jar dist/samreadlengthdistribution.jar src/test/resources/S*.bam

#ReadLength     S1      S2      S3      S4      S5
[-Inf/0[	0	0	0	0	0
[0/10[	0	0	0	0	0
[10/20[	0	0	0	0	0
[20/30[	0	0	0	0	0
[30/40[	0	0	0	0	0
[40/50[	0	0	0	0	0
[50/100[	1998	1998	1998	1998	1998
[100/150[	0	0	0	0	0
[150/200[	0	0	0	0	0
[200/250[	0	0	0	0	0
[250/300[	0	0	0	0	0
[300/350[	0	0	0	0	0
[350/400[	0	0	0	0	0
[400/450[	0	0	0	0	0
[450/500[	0	0	0	0	0
[500/1000[	0	0	0	0	0
[1000/Inf[	0	0	0	0	0
```


