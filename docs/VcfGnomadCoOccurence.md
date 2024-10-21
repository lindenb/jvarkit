# VcfGnomadCoOccurence

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Programmatic use of gnomad co-ocurence


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfgnomadcoocurence  [options] Files

Usage: vcfgnomadcoocurence [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
  * -B, --bed, --regions
      Bed intervals to scan
    --sleep
      Sleep 'x' seconds between each call
      Default: 1
    --sleep2
      Sleep 'x' seconds on error
      Default: 120
    --version
      print version and exit

```


## Keywords

 * vcf
 * gnomad



## Creation Date

20241018

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/coocurence/VcfGnomadCoOccurence.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/coocurence/VcfGnomadCoOccurence.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomadcoocurence** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 

## Example

```
java -jar dist/jvarkit.jarvcfgnomadcoocurence \
 src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz \
 --bed in.bed

(...)
    {
        "data": {
            "variant_cooccurrence": {
                "variant_ids": [
                    "1-905958-G-T",
                    "1-905962-C-A"
                ],
                "genotype_counts": [
                    124777,
                    1,
                    0,
                    8,
                    0,
                    0,
                    0,
                    0,
                    0
                ],
                "haplotype_counts": [
                    249563,
                    1,
                    8,
                    0
                ],
                "p_compound_heterozygous": 1,
                "populations": [
                    {
                        "id": "afr",
                        "genotype_counts": [
                            8042,
                            0,
                            0,
                            8,
                            0,
                            0,
                            0,
                            0,
                            0
                        ],
                        "haplotype_counts": [
                            16092,
                            0,
                            8,
                            0
                        ],
                        "p_compound_heterozygous": null
                    },
(...)

```

