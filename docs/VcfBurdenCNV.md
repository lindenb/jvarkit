# VcfBurdenCNV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Burden


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfburdencnv  [options] Files

Usage: vcfburdencnv [options] Files
  Options:
    --cases
      File or comma-separated list of control samples
    --controls
      File or comma-separated list of control samples
    --default-overlap
      default fraction overlap for option --known
      Default: 0.75
    --exclude-bed
      BED regions to exclude
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --include-bed
      BED regions to include
    --jexl
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatk.broadinstitute.org/hc/en-us/articles/360035891011 
    --known
      BED file containing the known frequent CNV. The 4th column must be a 
      number between 0 and 1, which well be the mutual overlap between the 
      known (e.g gnomad) variant and the user's variant. Default  value for 
      column 4 uses --default-overlap
    -o, --output
      Output file. Optional . Default: stdout
    --treshold
      do not display BED line if fisher > 'treshold'
      Default: 1.0E-5
    --vcfs
      path to vcfs. VCF must be indexed, one sample per vcf.
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * cnv



## Creation Date

20250404

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burdencnv/VcfBurdenCNV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burdencnv/VcfBurdenCNV.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdencnv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Still in beta. Do not use.



