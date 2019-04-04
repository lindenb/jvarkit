# SamShortInvertion

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan short inversions in SAM


## Usage

```
Usage: samshortinvert [options] Files
  Options:
    -B, --bed
      limit to that bed file
    -x, --extend
      extends interval by 'x' pb before merging.
      Default: 50
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --maxsize
      max size of inversion
      Default: 10000
    -o, --output
      Output file. Optional . Default: stdout
    -partition, --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -r, --rgn
      limit to that region CHROM:START-END
    --version
      print version and exit
    -s, -supporting
      Don't print the variant if INFO/DP <= 's'
      Default: 1

```


## Keywords

 * sam
 * bam
 * sv
 * inversion


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samshortinvert
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamShortInvertion.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamShortInvertion.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samshortinvert** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

finds the regions having some short inversions.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
```

