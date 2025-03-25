# RegenieMakeAnnot

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Create annotation files for regenie from a TSV input file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar regeniemakeannot  [options] Files

Usage: regeniemakeannot [options] Files
  Options:
    --gzip, -Z
      compress output files
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --masks
      mask file. TSV file. no header. 3 columns 
      prediction_name/score/comma-separated-mask_names. if undefined, will 
      produce one mask per prediction
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
  * -o, --output
      output dir.
    --prefix
      prefix for output files
      Default: chunk
    --reserve
      reserve 'n' output files of non-overlaping gene/target
      Default: 20
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --vcf
      keep variant id only if it's ID was found in this VCF file
    --version
      print version and exit
    -N
      max item per chunck. -1: no limit
      Default: 10000

```


## Keywords

 * vcf
 * regenie
 * burden



## Creation Date

20250311

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieMakeAnnot.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieMakeAnnot.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **regeniemakeannot** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


