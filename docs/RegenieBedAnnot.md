# RegenieBedAnnot

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Create annotation files for regenie using sliding annotations


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar regeniebedannot  [options] Files

Usage: regeniebedannot [options] Files
  Options:
    -A, --annotation
      value for annotation field
      Default: <empty string>
  * -B, --bed
      custom bed file chrom/start/end/name[/score]
    --chrom
      process only that chromosome
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --min-length
      slop each BED records in 5' and 3' so the minimal LENGTH is 'm'. 
      Multiple are comma separated
      Default: 0
    --noXY
      skip X/Y chromosome
      Default: false
    --version
      print version and exit
    -f
      comma separated of Allele frequencies , I will use the highest to 
      discard frequent variants.
      Default: 0.01
    -o
      Output file. Optional . Default: stdout

```


## Keywords

 * vcf
 * regenie
 * burden



## Creation Date

20250311

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieBedAnnot.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieBedAnnot.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **regeniebedannot** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





