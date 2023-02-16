# GtexRsToQTL

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

extract gtex eqtl data from a list of RS


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gtexrs2qtl  [options] Files

Usage: gtexrs2qtl [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gtex
 * rs
 * eqtl
 * sqtl



## Creation Date

20230215

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtex/GtexRsToQTL.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtex/GtexRsToQTL.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtexrs2qtl** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Motivation

export gtex eqtl data from a list of RS using the GTEX API

## Example

$ cat mylist.of.rs.txt | java -jar dist/jvarkit.jar gtexrs2qtl  | head | column -t

method            chromosome  datasetId  gencodeId           geneSymbol  geneSymbolUpper  nes        pValue       pos       snpId      tissueSiteDetailId                   variantId               phenotypeId
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.292217  1.86762e-15  35045523  rs6088690  Heart_Left_Ventricle                 chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101000.5   PROCR       PROCR            -0.341834  9.7875e-10   35045523  rs6088690  Heart_Left_Ventricle                 chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.360261  6.66882e-07  35045523  rs6088690  Brain_Putamen_basal_ganglia          chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101460.12  MAP1LC3A    MAP1LC3A         0.207036   3.12021e-10  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000078814.15  MYH7B       MYH7B            0.217765   3.34276e-17  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000100991.11  TRPC4AP     TRPC4AP          -0.174607  4.81976e-12  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.265171  2.8512e-15   35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101000.5   PROCR       PROCR            -0.181498  9.80399e-07  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000078814.15  MYH7B       MYH7B            0.211821   3.58767e-08  35045523  rs6088690  Esophagus_Gastroesophageal_Junction  chr20_35045523_A_G_b38  .


