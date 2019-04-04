# Biostar154220

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Cap BAM to a given coverage


## Usage

```
Usage: biostar154220 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -n, --depth
      number of reads
      Default: 20
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * coverage



## See also in Biostars

 * [https://www.biostars.org/p/154220](https://www.biostars.org/p/154220)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar154220
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar154220** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$ java -jar dist/sortsamrefname.jar --samoutputformat BAM input.bam  |\
  java -jar dist/biostar154220.jar  -n 20 --samoutputformat BAM |\
  samtools sort -T tmp -o output.bam -


$ samtools mpileup output.bam  | cut -f 4 | sort | uniq -c

  12692 0
 596893 1
  94956 10
  56715 11
  76947 12
  57912 13
  66585 14
  51961 15
  63184 16
  47360 17
  65189 18
  65014 19
 364524 2
 169064 20
  72078 3
 118288 4
  54802 5
  82555 6
  53175 7
  78474 8
  54052 9

```

