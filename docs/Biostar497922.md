# Biostar497922

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF into separate VCFs by SNP count. Deprecated use: VcfSplitNVariants


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar biostar497922  [options] Files

Usage: biostar497922 [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -C, --chrom
      by chromosomes when using '-n'.
      Default: false
    --count, -n
      number of variants per vcf
      Default: -1
    -d, --distance
      max distance beween consecutive variants (ignore if <=0) . A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --index
      start numering file index from 'x'
      Default: 1
    -D, --length
      max distance beween first and last variant (ignore if <=0) . A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    -m, --manifest
      output manifest to this file.
  * -o, --output
      Output directory
    --prefix
      file prefix
      Default: split
    --version
      print version and exit

```


## Keywords

 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/497922](https://www.biostars.org/p/497922)



## Creation Date

20210319

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar497922.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar497922.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar497922** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Deprecated

use **VcfSplitNVariants**

## Example

```
$ java -jar dist/biostar497922.jar -n 10 -o TMP src/test/resources/rotavirus_rf.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000001.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000002.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000003.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000004.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000005.vcf.gz
[INFO][Biostar497922]. Completed. N=45. That took:0 second


t$ find TMP/ -type f -name "*.vcf.gz" | sort | while read F; do echo -n "$F " && gunzip -c $F | grep -v "#" | wc -l  ; done
TMP/split.000001.vcf.gz 10
TMP/split.000002.vcf.gz 10
TMP/split.000003.vcf.gz 10
TMP/split.000004.vcf.gz 10
TMP/split.000005.vcf.gz 5

```

## Cited in

  * Yardeni G, Barfuss MHJ, Till W, Thornton MR, Crego CG, Lexer C, Leroy T, Paun O. The explosive radiation of the Neotropical Tillandsia subgenus Tillandsia (Bromeliaceae) has been accompanied by pervasive hybridization. Syst Biol. 2025 Jun 26:syaf039. doi: 10.1093/sysbio/syaf039. Epub ahead of print. PMID: 40569662.



