# AlmostSortedVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sort an 'almost' sorted VCF. Most variants should be sorted but a few  consecutive lines might have been switched by a caller.


## DEPRECATED

USE BCFTOOLS

## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar almostsortedvcf  [options] Files

Usage: almostsortedvcf [options] Files
  Options:
    -N, --buffersize
      max records in ram
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/AlmostSortedVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/AlmostSortedVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **almostsortedvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Warning

Barey used. Just use bcftools sort...

## Example

```bash
$  java -jar dist/almostsortedvcf.jar < in.vcf > out.vcf
```

