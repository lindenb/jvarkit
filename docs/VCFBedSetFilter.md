# VCFBedSetFilter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Set FILTER for VCF if intersects with BED.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfbedsetfilter  [options] Files

Usage: vcfbedsetfilter [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -e, --exclude, --blacklist
      Tribble or Tabix bed file containing the regions to be FILTERED. Must be 
      indexed with tribble or tabix, or use '--fast' to load in memory.
    -x, --extend
      Extend the variant coordinates per 'x' bases. A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 0
    -f, --filter
      FILTER name. If `--filter` is empty, FILTERED variant will be discarded.
      Default: VCFBED
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --fast, --memory
      Load the bed in memory: faster than tribble/tabix but memory consumming)
      Default: false
    --min-bed-fraction
      Min BED fraction overlap after extension. Only consider BED records if 
      variant overlap >= 'x' percent of bed length. A decimal number between 
      0.0 and 1.0. If the value ends with '%' it is interpretted as a 
      percentage eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. 
      e.g: '1/100' => '0.01'.
    --min-vc-fraction
      Min Variant fraction overlap after extension. Only consider BED records 
      if bed overlap >= 'x' percent of vc length. A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
    -M, --mutual-fraction
      Mutual comparator Two SV have are the same if they share a fraction 'x' 
      of their bases. For very small SV the fraction can be quite small while 
      for large SV the fraction should be close to 1. The Syntax is the 
      following : (<MAX_SIZE_INCLUSIVE>:<FRACTION as double or percent>;)+ . 
      For example if the SV as a size of 99bp, the fraction used with be 0.6 
      for '10:0.1;100:0.6;1000:0.9'. For the smallest size, a simple overlap 
      is a positive match.. Example:"10:0.5;100:75%;1000:80%;10000:90%".
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -i, --include, --whitelist
      Tribble or Tabix bed file containing the regions to be accepted. Regions 
      NOT overlapping those regions will be FILTERED. Must be indexed with 
      tribble or tabix, or use '--fast' to load in memory.

```


## Keywords

 * vcf
 * bed
 * filter



## See also in Biostars

 * [https://www.biostars.org/p/9465226](https://www.biostars.org/p/9465226)



## Creation Date

20150415

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbedsetfilter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Examples

```
$java -jar jvarkit.jar vcfbedsetfilter -f MYFILTER --include database.bed in.bed in.vcf
```



## Cited in

 * Megquier K, Turner-Maier J, Morrill K, Li X, Johnson J, Karlsson EK, et al. (2022) The genomic landscape of canine osteosarcoma cell lines reveals conserved structural complexity and pathway alterations. PLoS ONE 17(9): e0274383. https://doi.org/10.1371/journal.pone.0274383
 * Hannah Melhorn, Megan Sha, Manabu Kurihara, Katherine Megquier, Riley Aronson, Vicky Yang, Incidentally identified aortic dissection in a Great Pyrenees, Journal of Veterinary Cardiology, 2025, , ISSN 1760-2734, https://doi.org/10.1016/j.jvc.2025.09.007. (https://www.sciencedirect.com/science/article/pii/S1760273425001079)


