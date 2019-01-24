# VCFPolyX

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Number of repeated REF bases around POS.


## Usage

```
Usage: vcfpolyx [options] Files
  Options:
    -n, --filter
      if number of repeated bases is greater or equal to 'n' set a FILTER = 
      (tag) 
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
  * -R, --reference
      Indexed Genome Reference. It can be a the path to fasta file that must 
      be indexed with samtools faidx and with picard CreateSequenceDictionary. 
      It can also be a BioDAS dsn url like 
      `http://genome.cse.ucsc.edu/cgi-bin/das/hg19/` . BioDAS references are 
      slower, but allow to work without a local reference file.
    --skip-filtered
      Don't spend some time to calculate the tag if the variant is FILTERed
      Default: false
    -t, --tag
      Tag used in INFO and FILTER columns.
      Default: POLYX
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * repeat


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfpolyx
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyX.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyX.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyXTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyXTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpolyx** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java  -jar dist/vcfpolyx.jar -R reference.fa input.vcf
(...)
2   1133956 .   A   G   2468.84 .   POLYX=23
2   1133956 .   A   AG  3604.25 .   POLYX=23
2   2981671 .   T   G   47.18   .   POLYX=24
(...)
```

