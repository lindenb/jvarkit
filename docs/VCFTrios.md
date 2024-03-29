# VCFTrios

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find mendelian incompatibilitie / denovo variants in a VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcftrio  [options] Files

Usage: vcftrio [options] Files
  Options:
    -A, --attribute
      INFO Attribute name containing the name of the affected samples.
      Default: MENDEL
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --discard
      Discard the variant if there is NO mendelian violation.
      Default: false
    -fi, --filter-in
      FILTER name if there is ANY mendelian violation.
    -fo, --filter-out, --filter-no-denovo
      FILTER name if there is NO mendelian violation.
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -gtf, --gtfilter
      GENOTYPE FILTER name. Create a filter in the GENOTYPE column when there 
      is NO mendelian violation
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hr, --hom-ref
      [20180705] treat NO_CALL genotypes as HOM_REF (when individual 
      VCF/Sample have been merged).
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -p, --ped, --pedigree
      Pedigree file. A pedigree file. tab delimited. Columns: 
      family,id,father,mother, 
      sex:(0|.|undefined|unknown:unknown;1|male|M:male;2|female|F:female), 
      phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --version
      print version and exit

```


## Keywords

 * vcf
 * mendelian
 * pedigree
 * denovo



## Creation Date

20130705

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcftrios/VCFTrios.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcftrios/VCFTrios.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcftrios/VCFTriosTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcftrios/VCFTriosTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcftrio** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Using GATK VariantAnnotator ? 

GATK VariantAnnotator doesn't work if GQ is low or if there is no GQ.

## Example

a pedigree file:

```
$  cat pedigree.txt 

A	SAMPLE_P	0	0	0
A	SAMPLE_M	0	0	0
A	SAMPLE_E	SAMPLE_P	SAMPLE_M	0
```


find mendelian incompatibilities:

```
$  gunzip -c input.vcf.gz |\
   java -jar dist/vcftrio.jar -p pedigree.txt | grep -E '(#CHROM|MENDEL=SAMPLE_E)' |\
   verticalize 

(...)
>>> 23
$1	#CHROM	X
$2	POS	0573
$3	ID	rs358
$4	REF	G
$5	ALT	A
$6	QUAL	85.60
$7	FILTER	PASS
$8	INFO	MENDEL=SAMPLE_E
$9	FORMAT	GT:DP:DP4:GP:GQ:PL
$10	SAMPLE_E	0/1:11:6,0,5,0:97,0,122:97:96,0,118
$11	SAMPLE_M	1/1:5:0,0,5,0:134,19,0:19:120,15,0
$12	SAMPLE_P	1/1:6:0,0,6,0:136,22,0:22:121,18,0
<<< 23
(...)
>>> 59
$1	#CHROM	Y
$2	POS	19
$3	ID	rs5678
$4	REF	CA
$5	ALT	C,CAA
$6	QUAL	31.86
$7	FILTER	PASS
$8	INFO	MENDEL=SAMPLE_E
$9	FORMAT	GT:DP:DP4:GP:GQ
$10	SAMPLE_E	2/2:80:3,0,43,34:.,.,108,.,203,0:99
$11	SAMPLE_M	.
$12	SAMPLE_P	1/1:53:0,0,27,26:81,99,0,.,.,.:81
<<< 59

```
## History

  * [20180907] moved to the new DeNovoDetector
  * [20180704] changing the arguments that are not really clear.

 

