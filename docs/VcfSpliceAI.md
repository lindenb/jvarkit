# VcfSpliceAI

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate VCF with local spiceai vcf


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfspliceai  [options] Files

Usage: vcfspliceai [options] Files
  Options:
  * --vcf, --annotation, --spliceai
      SpliceAI VCF.vcf.gz indexed with tabix
      Default: []
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --buffer
      When we're looking for variants in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 1000
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --tag
      INFO tag
      Default: SpliceAI
    --version
      print version and exit

```


## Keywords

 * vcf
 * splice
 * splicing
 * spliceai



## Creation Date

20201107

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/spliceai/VcfSpliceAI.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/spliceai/VcfSpliceAI.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfspliceai** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Example

```
java -jar dist/jvarkit.jar vcfspliceai --annot /path/to/spliceai_scores.masked.indel.hg38.vcf.gz  src/test/resources/test_vcf01.vcf 

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
(...)
1	866893	.	T	C	431	PASS	AA=t;AC=7;AF=0.7;AN=10;SpliceAI=SAMD11|0.00|0.00|0.00|0.00|-13|24|-13|-48
1	870317	.	G	A	12	PASS	AC=11;AF=0.917;AN=12;SpliceAI=SAMD11|0.00|0.00|0.00|0.00|2|17|16|-12
1	875770	.	A	G	338	PASS	AA=a;AC=8;AF=0.8;AN=10;SpliceAI=SAMD11|0.00|0.00|0.01|0.00|-1|-45|-50|-46
1	903245	.	A	G	199	PASS	AA=a;AC=6;AF=0.6;AN=10;SpliceAI=PLEKHN1|0.00|0.00|0.00|0.00|48|-37|-22|1
1	905130	.	ATG	A	487	PASS	AC=3;AF=0.5;AN=6;CIGAR=1M2D;IDREP=1;REFREP=2;RU=TG;SpliceAI=PLEKHN1|0.00|0.00|0.00|0.00|-43|21|-33|-37
1	909238	.	G	C	229	PASS	AA=C;AC=8;AF=0.667;AN=12;SpliceAI=PLEKHN1|0.00|0.01|0.00|0.00|-43|-50|39|-7
1	912049	.	T	C	400	PASS	AA=T;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.01|0.01|0.00|-28|-14|-27|-23
1	913889	.	G	A	372	PASS	AA=G;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.01|0.00|0.00|-46|9|2|-45
1	914333	.	C	G	556	PASS	AA=G;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|-3|27|-3|-38
1	914852	.	G	C	525	PASS	AA=C;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|22|-22|48|49
1	914940	.	T	C	488	PASS	AA=C;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|28|-30|-39|3
(...)
```


