# SamStructVarNeg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find Structural Variation by Negative Comparaison


## Usage

```
Usage: svneg [options] Files
  Options:
    -B, --bed
      Optional BED file. SV should overlap this bed.
    --chrom-regex
      Only consider the chromosomes matching the following regular expression.
      Default: (chr)?([1-9][0-9]*|X|Y)
  * -b, --controls, --bams
      Control Bams. One path per lines
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
    --min-controls
      Min number controls bam matching the event
      Default: 1
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
$ ./gradlew svneg
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190413

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNeg.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNeg.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNegTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNegTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **svneg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



