# SamToJson

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert a SAM input to JSON


## Usage

```
Usage: sam2json [options] Files
  Options:
    -atts, --atts
      do not print attributes
      Default: false
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -cigar, --cigar
      expand cigar
      Default: false
    -flag, --flag
      expand SAm Flags
      Default: false
    -H, --header
      don't print SAM HEADER
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -name, --name
      do not print read name
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * json


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sam2json
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210402

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToJson.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToJson.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam2json** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



