# VcfGatkEval

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Eval/Plot gatk INFO tags for filtering


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfgatkeval  [options] Files

Usage: vcfgatkeval [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -I, --input-type
      input type. vcf: vcf file or stdin. table: stdin or one or more output 
      of *.output.table.txt
      Default: vcf
      Possible Values: [vcf, table]
    -o, --output
      filename prefix
      Default: gatk.eval
    -p, --percentile
      GATK Filters should be applied to this percentile: f < x < (1.0 -f)
      Default: 0.025
    --version
      print version and exit
    --depth, --with-depth
      include INFO/DP
      Default: false

```


## Keywords

 * vcf
 * gatk



## Creation Date

20230424

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgatkeval/VcfGatkEval.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgatkeval/VcfGatkEval.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgatkeval** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
# Run on one vcf

```
$ bcftools view in.bcf | java -jar dist/jvarkit.jar  vcfgatkeval -o "out1"
```

list output

```
$ ls out1.*
out1.output.filters.txt
out1.output.R 
out1.output.table.txt
```

plot barplots:

```
R --vanilla --no-save < out1.output.R 
```


filters for gatk that can be used using gatk `--arguments_file`

```
cat  out1.output.filters.txt
```

```
-filter
"vc.isSNP() && FS > 21.0"
--filter-name
FS_HIGH_SNP
-filter
"vc.isSNP() && MQ < 60.0"
--filter-name
MQ_LOW_SNP
-filter
"vc.isSNP() && MQRankSum < 0.0"
--filter-name
MQRankSum_LOW_SNP
-filter
"vc.isSNP() && MQRankSum > 0.0"
--filter-name
MQRankSum_HIGH_SNP
-filter
"vc.isSNP() && QD < 1.0"
--filter-name
QD_LOW_SNP
-filter
"vc.isSNP() && ReadPosRankSum < -2.2"
--filter-name
ReadPosRankSum_LOW_SNP
-filter
"vc.isSNP() && ReadPosRankSum > 2.4"
--filter-name
ReadPosRankSum_HIGH_SNP
-filter
"vc.isSNP() && SOR > 3.5"
--filter-name
SOR_HIGH_SNP
```

use with gatk variantFilteration:

```
gatk VariantFiltration -V in.vcf.gz -R reference.fasta -O out.vcf.gz --arguments_file out1.output.filters.txt

```

# Parallelisation

run in parallel
```
$ bcftools view in.bcf chr1 | java -jar dist/jvarkit.jar  vcfgatkeval -o "out1"
$ bcftools view in.bcf chr2 | java -jar dist/jvarkit.jar  vcfgatkeval -o "out2"
```

and then concat:

```
cat out1.output.table.txt out2.output.table.txt  |  java -jar dist/jvarkit.jar  vcfgatkeval --input-type table  -o "out3"
```



