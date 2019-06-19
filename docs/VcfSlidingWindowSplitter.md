# VcfSlidingWindowSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF by sliding window


## Usage

```
Usage: vcfwindowsplitter [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -l, --list
      list all available extractors
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -n, --min-variant
      Minimum number of variants to write a vcf. don't write if num(variant) < 
      'x' 
      Default: 1
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -s, -S, --window-shift
      Sliding window shift. A distance specified as a positive integer.Comma 
      are removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 500000
    -w, -W, --window-size
      Sliding window size. A distance specified as a positive integer.Comma 
      are removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000000

```


## Keywords

 * vcf
 * sliding
 * window


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfwindowsplitter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190619

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitter.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfwindowsplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
$ java -jar dist/vcfwindowsplitter.jar -w 1000 -s 500 -o jeter.zip -m jeter.manifest src/test/resources/rotavirus_rf.vcf.gz 


```

