# VCFNearest

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

find nearest feature near a variant


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfnearest  [options] Files

Usage: vcfnearest [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
  * -B, --bed
      Bed path.
  * -C, --columns
      comma separated names of the columns in the bed file starting from the 
      4th column. Empty columns will be skipped. e.g: 'NAME,TYPE,,SCORE' 
      expects a bed with 7 columns. The 6th will be skipped. Two columns are 
      added: the side (0: overlap, -1 bed is before variant , 1 bed is after 
      variant ) and the distance
      Default: <empty string>
    -d, --distance
      max distance from variant to feature. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -n, --max-features
      do not print more than 'x' features; disable if x <=0
      Default: -1
    -o, --out
      Output file. Optional . Default: stdout
    -t, --tag
      Tag name
      Default: NEAREST
    --version
      print version and exit

```


## Keywords

 * vcf
 * vep
 * annotation



## Creation Date

20250523

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfnearest/VCFNearest.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfnearest/VCFNearest.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfnearest/VCFNearestTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfnearest/VCFNearestTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfnearest** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


