# FindGVCFsBlocks

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find common blocks of calleable regions from a set of gvcfs


## Usage

```
Usage: findgvcfsblocks [options] Files
  Options:
    --min-size, --block-size
      min block size. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 0
    -c, --chrom, --chromosome, --contig
      limit to that contig
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -T
      temporary directory

```


## Keywords

 * gvcf
 * gatk
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew findgvcfsblocks
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210806

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gvcf/FindGVCFsBlocks.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gvcf/FindGVCFsBlocks.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findgvcfsblocks** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

find regions for running GATK CombineGVCFs in parallel.

## Input

input is a set of path to the indexed g.vcf files
or it's a file with the '.list' suffix containing the path to the g.vcf files

g.vcf files must be indexed if option `-c` is used.

## Output

output is a BED file containing the calleable GVCFs blocks.

## Example

```
$ java -jar dist/findgvcfsblocks.jar --chrom RF11 S1.g.vcf.gz S2.g.vcf.gz S3.g.vcf.gz 
RF11	0	5
RF11	5	12
RF11	12	15
RF11	15	18
RF11	18	20
RF11	20	21
RF11	21	27
RF11	27	28
RF11	28	30
RF11	30	46
(...)
```
