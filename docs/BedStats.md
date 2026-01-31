# BedStats

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

statistics about one or more file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedstats  [options] Files

Usage: bedstats [options] Files
  Options:
    --disable-normalize-contig
      do not normalize contig name '1'->'chr1'
      Default: false
    -f, --fraction
      for overlap between to BED file. Two BED record overlap if they  both 
      share 'x' fraction of overlap compared to their lengthA decimal number 
      between 0.0 and 1.0. If the value ends with '%' it is interpretted as a 
      percentage eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. 
      e.g: '1/100' => '0.01'.
      Default: 1.0E-7
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mqc_id
      mqc id
      Default: bedstats_
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * stats
 * multiqc



## Creation Date

20260130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedstats/BedStats.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedstats/BedStats.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bedstats/BedStatsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bedstats/BedStatsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedstats** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


input is bed on standard input or it's a set of interval files (.bed, .interval_list, .gtf, etc... )

output is text file containing multiple chuncks for multiqc


## Example

```
$ java -jar dist/bedstats.jar < in1.bed > out.txt
$ java -jar dist/bedstats.jar in1.bed in2.bed in3.bed > out.txt
```



