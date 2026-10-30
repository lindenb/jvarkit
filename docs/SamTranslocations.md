# SamTranslocations

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Explore balanced translocations between two chromosomes using discordant paired-end reads.


## Usage

```
Usage: java -jar dist/samtransloc.jar  [options] Files
Usage: samtransloc [options] Files
  Options:
    -B, --bed
      Optional BED file. SV should overlap this bed.
    --chrom-regex
      Only consider the chromosomes matching the following regular expression.
      Default: (chr)?[0-9XY]+
    -d, --distance
      Max distance between two read to test if they both end at the same ~ 
      position. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality.
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --version
      print version and exit
    -exclude-bed
      skip reads pverlapping that bed file

```


## Keywords

 * sam
 * bam
 * sv
 * translocation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samtransloc
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20171130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samtransloc/SamTranslocations.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samtransloc/SamTranslocations.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samtransloc/SamTranslocationsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samtransloc/SamTranslocationsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samtransloc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a set of bam files or one file with the suffix '.list' containing the path to the bams.

## Example

```
$ java -jar dist/jvarkit.jar samtransloc src/test/resources/HG02260.transloc.chr9.14.bam  > out.json
```


