# SamTranslocations

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Explore balanced translocations between two chromosomes using discordant paired-end reads.


## Usage

```
Usage: samtranslocations [options] Files
  Options:
    -B, --bed
      Optional BED file. SV should overlap this bed.
    --chrom-regex
      Only consider the chromosomes matching the following regular expression.
      Default: (chr)?([1-9][0-9]*|X|Y)
    -d, --distance
      Max distance between two read to test if they both end at the same ~ 
      position. A distance specified as a positive integer.Comma are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --low-complexity
      zone of low complexity (=many discordant reads). In those region the 
      software will 'freeze'. If the number of buffered reads * this number 
      per sample then clear the buffer
      Default: 2000
    --mapq
      min mapping quality.
      Default: 0
    --max-sa
      ignore reads having more that SA:Z supplementary alignments.
      Default: 4
    -m, --min
      Min number of events to validate the translocation
      Default: 3
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * sv
 * translocation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samtranslocations
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamTranslocations.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamTranslocations.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamTranslocationsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamTranslocationsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samtranslocations** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a set of bam files or one file with the suffix '.list' containing the path to the bams.

## Example

```
$ java -jar dist/samtranslocations.jar src/test/resources/HG02260.transloc.chr9.14.bam | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG02260
9	137230996	9:137230996:14:79839048	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=14;DP=18;POS2=79839048;STDDEV_POS1=120;STDDEV_POS2=187;SVTYPE=BND	GT:DP:SR	0/1:18:5,13,13,5
14	79839131	14:79839131:9:137230969	N	<TRANSLOC>	17	.	AC=1;AF=0.500;AN=2;CHROM2=9;DP=17;POS2=137230969;STDDEV_POS1=153;STDDEV_POS2=153;SVTYPE=BND	GT:DP:SR	0/1:17:12,5,5,12
```

## History

* 2018-09-18 :  rewriting
* 2017-12-13 :  refactoring for balanced translocation.

