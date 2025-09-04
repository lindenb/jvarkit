# VcfISec

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Only prints variants that are contained/not contained into another VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfisec  [options] Files

Usage: vcfisec [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --vcf, --vcfs, --file-list
      External indexed VCFs. a file with the suffix '.list' is interpreted as 
      a file containing the path to the indexed VCFs
      Default: []
    --filter
      FILTER name for SOFT filtering
      Default: <empty string>
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --inverse
      Print variants that are not part of the VCF-database.
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
    --version
      print version and exit

```


## Keywords

 * vcf
 * compare



## See also in Biostars

 * [https://www.biostars.org/p/287815](https://www.biostars.org/p/287815)



## Creation Date

20140204

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfisec/VcfISec.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfisec/VcfISec.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfisec/VcfISecTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfisec/VcfISecTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfisec** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


* For Structural Variants (identified with INFO/SVTYPE) we just compare the fraction of common overlap.
* For SNVs : spanning deletions (`*`) are not considered when comparing alleles

# Example

```
find DIR -name "*.vcf.gz" > paths.list
java -jar jvarkit.jar vcfisec --vcfs paths.list in.vcf
```

# See also:

 * bcftools isec


