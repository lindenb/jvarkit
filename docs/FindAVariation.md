# FindAVariation

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Finds a specific mutation in a list of VCF files


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar findavariation  [options] Files

Usage: findavariation [options] Files
  Options:
    --enable-non-indexed, --no-index
      Enable Scanning of non-indexed VCF/BCF files
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -homref, --homref
      Hide HOM_REF genotypes
      Default: false
    -nocall, --nocall
      Hide NO_CALL genotypes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --snv-matcher
      How to test if two SNVs are the same. Spanning deletion will be ignored.
      Default: chrom_pos_ref_any_alt
      Possible Values: [id, overlap, chrom_pos, chrom_pos_ref, chrom_pos_ref_any_alt, chrom_pos_ref_all_alt]
    --sv-fraction
      How to match two SV. Two SV have are the same if they share a fraction 
      'x' of their bases. For very small SV the fraction can be quite small 
      while for large SV the fraction should be close to 1. The Syntax is the 
      following : (<MAX_SIZE_INCLUSIVE>:<FRACTION as double or percent>;)+ . 
      For example if the SV as a size of 99bp, the fraction used with be 0.6 
      for '10:0.1;100:0.6;1000:0.9'. For the smallest size, a simple overlap 
      is a positive match.
      Default: 10:0.5;100:0.75;1000:0.8;10000:0.9
  * -V, --vcf, --variants
      A VCF file containing the variants
    --version
      print version and exit

```


## Keywords

 * vcf
 * variation
 * search
 * find
 * bcf



## Creation Date

20140623

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/findavariation/FindAVariation.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/findavariation/FindAVariation.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/findavariation/FindAVariationTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/findavariation/FindAVariationTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findavariation** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

```
$ find ./ -name "*.vcf" -o -name "*.vcf.gz" |\
   java -jar dist/findavariation.jar -V search.vcf


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

 

