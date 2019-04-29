# MakeMiniBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates an archive of small bams with only a few regions.


## Usage

```
Usage: mkminibam [options] Files
  Options:
    -B, --bed
      Use the intervals from this BED file.
    --bnd
      [20190427]When reading VCF file, don't get the mate position for the 
      structural BND variants.
      Default: false
    -b, --bounds, --edge
      [20190427] If `b` is greater than 0 and the genomic object has a length 
      greater than `b` then consider the edges of the object as two positions. 
      the idea is to just save the boundaries of a large deletion. A distance 
      specified as a positive integer.Comma are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: -1
    -C, --comment
      [20190427]Add a file '*.md' with this comment.
      Default: <empty string>
    -x, --extend
      Extend the positions by 'x' bases. A distance specified as a positive 
      integer.Comma are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 5000
    --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: Accept All/ Filter out nothing
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
      File prefix in the archive. Special value 'now' will be replace by the 
      current date
      Default: miniBam.
    -R, --reference
      Optional Reference file for CRAM files. Indexed fasta Reference file. 
      This file must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary 
    -T, --tmp
      Tmp working directory
      Default: /tmp
    -V, --variant, --vcf
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

# Input

 input is a set of bam files or a file with the suffix '.list' containing one path to a bam per line.
 
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
 
 
