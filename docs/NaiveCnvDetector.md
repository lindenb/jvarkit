# NaiveCnvDetector

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

experimental CNV detection for multiple samples.


## Usage

```
Usage: naivecnvdetector [options] Files
  Options:
    -c, --config
      config file. Tab delimited. 
      Sample-name(tab)mean-depth(tab)integer[affected=1,non-affected=0]. If 
      this file is not specified , all samples are considered unaffected 
      (discovery mode).
    -E, -del, --del, --deletion
      Deletion Treshold. Which fraction of the median depth is considered as 
      aa deletion. Must be <1.0
      Default: 0.5
    --disable-consecutive
      Disable 'consecutive' positions criteria. Default: dump buffer of 
      position if the current samtools-depth line is not the very next 
      expected position.
      Default: false
    -U, -dup, --dup, --duplication
      Duplication Treshold. Which fraction of the median depth is considered 
      as a duplication. Must be >1.0
      Default: 1.5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -md, --min-dp
      At least one 'unaffected' sample must have a normalized-depth greater 
      than this value.
      Default: 20
    --no-both
      There cannot be a DEL and a DUP at the same place.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -stddevu, --stddev-unaffected
      Maximum standard deviation of depth for unaffected samples. Ignored if 
      negative or if not any affected samples is defined.
      Default: 10.0
    --version
      print version and exit
    --weirdDepth
      Treat normalized depth greater than this value as 'weird' and discard 
      the sliding windows at this place.
      Default: 500
    -R, -reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -s
      window shift
      Default: 500
    -t
      DEL must be < median-depth-stdev and DUP must be > median-depth+stdev
      Default: false
    -w
      window size
      Default: 1000

```


## Keywords

 * cnv
 * bam
 * sam
 * wig
 * bigwig
 * bigbed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew naivecnvdetector
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/NaiveCnvDetector.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/NaiveCnvDetector.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/NaiveCnvDetectorTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/NaiveCnvDetectorTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **naivecnvdetector** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Input is either:

  * one fileof samtools depth output. All 'N' samples in one file.
  * 'N' files (samtools depth output AND/OR bigwig/bigbed [experimental not tested] ). One samples in per file. REF dictionary is required. List of file can be specified if input ends with '.list' 

bigbed and bigwig have not been tested; Bigbed shouldn't have overlapping regions...

## Example

```
samtools depth -r '1:1234-567' *.bam |\
	java -jar dist/naivecnvdetector.jar  > out.tsv
```


