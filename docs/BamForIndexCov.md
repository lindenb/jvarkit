# BamForIndexCov

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

prepare BAM/CRAM from indexcov.


## Usage

```
Usage: java -jar dist/bam4indexcov.jar  [options] Files
Usage: bam4indexcov [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --include-chromosomes
      regex of chromosomes to only include
      Default: (chr)?[0-9XY]+
    -m, --manifest
      manifest file
    -Q, --mapq
      min mapping quality.
      Default: 10
    --md5
      generate md5 file
      Default: false
  * -o, --output
      Output directory
    -p, --prefix
      File prefix
      Default: bam4indexcov.
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * cnv
 * duplication
 * deletion
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam4indexcov
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20220506

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/BamForIndexCov.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/BamForIndexCov.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam4indexcov** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## example

```
java -jar dist/bam4indexcov.jar -R src/test/resources/rotavirus_rf.fa -o  . src/test/resources/S*.bam -i 'RF1.*'
