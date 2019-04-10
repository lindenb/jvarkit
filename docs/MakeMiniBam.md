# MakeMiniBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates an archive of small bams with only a few region


## Usage

```
Usage: mkminibam [options] Files
  Options:
    -B, --bed
      Use the intervals from this BED file.
    -x, --extend
      Extend the positions by 'x' bases. A distance specified as a positive 
      integer.Comma are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 5000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    -p, --pos
      Add this position 'chrom:pos'
      Default: []
    --prefix
      File prefix in the archive
      Default: miniBam.
    -T, --tmp
      Tmp working directory
      Default: /tmp
    -V, --variant
      Use the intervals from this VCF file.
    --version
      print version and exit

```


## Keywords

 * bam
 * sam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew mkminibam
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190410

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/minibam/MakeMiniBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/minibam/MakeMiniBam.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/minibam/MakeMiniBamTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/minibam/MakeMiniBamTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **mkminibam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
# Motivation

Bams are too bigs and my users often ask to visualize a small region of a set of bam
 
# Example
 
```
$  find src/test/resources/ -name "S*.bam" > bams.list
$   java -jar dist/mkminibam.jar -p "RF01:100" -o out.zip bams.list 
[INFO][MakeMiniBam]src/test/resources/S5.bam
[INFO][MakeMiniBam]src/test/resources/S2.bam
[INFO][MakeMiniBam]src/test/resources/S4.bam
[INFO][MakeMiniBam]src/test/resources/S3.bam
[INFO][MakeMiniBam]src/test/resources/S1.bam

$ unzip -t out.zip 

Archive:  out.zip
    testing: miniBam.S5.bam           OK
    testing: miniBam.S5.bai           OK
    testing: miniBam.S2.bam           OK
    testing: miniBam.S2.bai           OK
    testing: miniBam.S4.bam           OK
    testing: miniBam.S4.bai           OK
    testing: miniBam.S3.bam           OK
    testing: miniBam.S3.bai           OK
    testing: miniBam.S1.bam           OK
    testing: miniBam.S1.bai           OK
No errors detected in compressed data of out.zip.
```
 
 
