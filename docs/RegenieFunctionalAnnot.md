# RegenieFunctionalAnnot

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Create annotation files for regenie using snpEff annotations


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar regeniefunctionalannot  [options] Files

Usage: regeniefunctionalannot [options] Files
  Options:
  * -A, --annotations
      seq_ontology <-> score file. TSV file. no header. at least 2 columns 
      prediction_name/score 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieFunctionalAnnot.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/regenie/RegenieFunctionalAnnot.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **regeniefunctionalannot** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

snpeff.cadd.in.vcf.gz |\
	java -jar dist/jvarkit.jar regeniefunctionalannot -A anot.tsv |\
	java -jar dist/jvarkit.jar regeniemakeannot -o OUT



