# CommBams

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Equivalent of unix 'comm' for bams sorted on queryname


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar commbams  [options] Files

Usage: commbams [options] Files
  Options:
    -delim, --delimiter
      Output delimiter
      Default: 	
    -empty, --empty
      Empty content symbol
      Default: .
    -f, --format
      What should I print ? (only the read name ? etc...)
      Default: name
      Possible Values: [name, but_metadata, all]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -1, --hide1
      suppress read unique to file 1
      Default: false
    -2, --hide2
      suppress read unique to file 2
      Default: false
    -3, --hide3
      suppress reads present in both files
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -sortmethod, --sortmethod
      [20171110]Method used to sort the read on query name. (samtools != 
      picard) see https://github.com/samtools/hts-specs/issues/5
      Default: picard
      Possible Values: [samtools, picard]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * comm
 * compare



## See also in Biostars

 * [https://www.biostars.org/p/9493549](https://www.biostars.org/p/9493549)



## Creation Date

20170420

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBams.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBams.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBamsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBamsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **commbams** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


