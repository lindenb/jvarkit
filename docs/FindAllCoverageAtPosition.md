# FindAllCoverageAtPosition


## Usage

```
Usage: findallcoverageatposition [options] Files
  Options:
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    --groupby
      Group Reads by
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    -o, --out
      Output file. Optional . Default: stdout
    -f, --posfile
      File containing positions. if file suffix is '.bed': all positions in 
      the range will be scanned.
    -p, --position
      -p chrom:pos . Multiple separated by space. Add this chrom/position. 
      Required 
      Default: <empty string>
    --version
      print version and exit

```


## Description

Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ?


## Keywords

 * bam
 * coverage
 * search
 * depth


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make findallcoverageatposition
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAllCoverageAtPosition.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAllCoverageAtPosition.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findallcoverageatposition** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

```
$ find ./  -type f -name "*.bam" |\
   java -jar dist/findallcoverageatposition.jar -p "chr2:1234" 


#File   CHROM   POS SAMPLE  DEPTH   M   I   D   N   S   H   P   EQ  X
/path/to/Sample1.bam    2   1234    SAMPLE1 10  10  0   1   0   0   0   0   0   0
/path/to/Sample2.bam    2   1234    SAMPLE2 10  0   0   0   1   0   0   0   5   5
/path/to/Sample3.bam    2   1234    SAMPLE3 10  10  0   0   0   0   0   0   0   0
```

## See also

 * https://twitter.com/pjacock/status/538300664334798848
 * https://twitter.com/yokofakun/status/538300434109456385
 * https://twitter.com/pjacock/status/538299549455233024
 * FindAVariation

## History

 * 2017: moved to jcommander


