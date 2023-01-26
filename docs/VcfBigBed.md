# VcfBigBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with values from a bigbed file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfbigbed  [options] Files

Usage: vcfbigbed [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
  * -B, --bed, --bigbed
      Path to the bigbed file.
    --bufferSize
      When we're looking for bed in a lare bigbed file, load the bed items in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 10000
    -e, --expr, --jexl, --format
      A JEXL Expression returning a string JEXL stands for Java EXpression 
      Language.  See 
      https://commons.apache.org/proper/commons-jexl/reference/syntax.html . 
      The variable 'bed' is the current observed Bed Line. It implements 
      java.util.List<String> and Locatable. The variable 'ctx' or 'variant' is 
      the current observed variant. The variable 'line' is the original bed 
      line 
      Default: bed.get(0)+":"+bed.get(1)+"-"+bed.get(2)
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -mofb, --min-overlap-bed-fraction
      Minimum overlap required as a fraction of BIGBED record. A decimal 
      number between 0.0 and 1.0. If the value ends with '%' it is 
      interpretted as a percentage eg. '1%' => '0.01'. A slash '/' is 
      interpretted as a ratio. e.g: '1/100' => '0.01'.
    -mofv, --min-overlap-vcf-fraction
      Minimum overlap required as a fraction of VCF record. A decimal number 
      between 0.0 and 1.0. If the value ends with '%' it is interpretted as a 
      percentage eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. 
      e.g: '1/100' => '0.01'.
    -o, --out
      Output file. Optional . Default: stdout
    -T, --tag, -tag
      Name of the INFO tag. default: name of the bigbed
    --version
      print version and exit

```


## Keywords

 * vcf
 * wig
 * wiggle
 * bigbed
 * bed



## Creation Date

20220107

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VcfBigBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VcfBigBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbigbed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ java -jar dist/vcfbigbed.jar \
	-B "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/JASPAR2022_hg19.bb" \
	--format  'bed.get(1)+"-"+bed.get(2)+":"+bed.get(6)'  \
	input.vcf

(...)
##INFO=<ID=JASPAR2022_hg19,Number=.,Type=String,Description="Values from bigbed file: http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/JASPAR2022_hg19.bb format bed.get(1)+\"-\"+bed.get(2)+\":\"+bed.get(6)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr22	41697508	.	A	C	48.67	.	AC=2;AN=10;BQB=0.572843;DP=36;DP4=19,7,3,5;HOB=0.32;ICB=0.425;JASPAR2022_hg19=41697498-41697509:ZBTB12,41697506-41697515:NKX2-8,41697507-41697524:Spi1,41697498-41697509:Stat5a::Stat5b,41697503-41697509:Foxn1,41697507-41697519:ZNF263,41697493-41697509:BCL6,41697503-41697517:TFAP2C,41697502-41697510:NR2C2;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.658863;SGB=10.3229;VDB=0.693968	GT:PL	0/0:0,9,47	0/0:0,18,73
```



