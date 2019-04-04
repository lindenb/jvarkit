# CytobandToSvg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates a svg karyotype .


## Usage

```
Usage: cytoband2svg [options] Files
  Options:
    -C, --cytobands
      Cytoband URI. tab delimited file. no header. 
      chrom/chromStart/chromEnd/name/gieStain. '-'=stdin
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
    -H, --height
      Image height
      Default: 700
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -t, --title
      title
      Default: <empty string>
    --version
      print version and exit
    -W, -width, --width
      Image width
      Default: 1000

```


## Keywords

 * karyotype
 * svg
 * ideogram


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew cytoband2svg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/CytobandToSvg.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/CytobandToSvg.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/CytobandToSvgTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/CytobandToSvgTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **cytoband2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a tab delimited file similar to http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

no header

columns are

* chrom
* chromStart
* chromEnd
* name
* gieStain

Chromosomes will be ordered occording to the first time they're seen in the file, so you'd better sort them.

```
$ wget -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" |\
	gunzip -c |\
	sort -t $'\t' -k1,1V -k2,2n > cytoBand.txt
```

