# VcfAlleleBalance

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Insert missing allele balance annotation using FORMAT:AD


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfallelebalance  [options] Files

Usage: vcfallelebalance [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -f, --filtered
      ignore FILTER-ed **GENOTYPES**
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -p, -ped, --pedigree, --ped
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0|.|undefined|unknown:unknown;1|male|M:male;2|female|F:female), 
      phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -s, --snp
      consider only snps
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * allele-balance
 * depth



## Creation Date

20180829

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/allelebalance/VcfAlleleBalance.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/allelebalance/VcfAlleleBalance.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/allelebalance/VcfAlleleBalanceTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/allelebalance/VcfAlleleBalanceTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfallelebalance** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

August 2018, for @MKarakachoff

part of the code was inspired from GATK public code :  https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/annotator/AlleleBalance.java

## Example

```
$ java -jar dist/vcfallelebalance.jar  src/test/resources/test_vcf01.vcf
$ java -jar dist/vcfallelebalance.jar -p src/test/resources/test_vcf01.ped src/test/resources/test_vcf01.vcf
```


