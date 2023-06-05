# VcfStats

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Produce VCF statitics


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfstats  [options] Files

Usage: vcfstats [options] Files
  Options:
    -exclude, --exclude
      name of modules to be excluded
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --list
      list available modules and exit
    -o, --output
      Output file. Optional . Default: stdout
    --prefix
      file prefix
      Default: <empty string>
    --version
      print version and exit
    -D
      other parameters.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * vcf
 * stats
 * R



## Creation Date

20131212

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstats** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/jvarkit.jar vcfstats src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz |  R --no-save 
```

