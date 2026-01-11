# VcfGnomadSV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek annotations from gnomad structural variants


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfgnomadsv  [options] Files

Usage: vcfgnomadsv [options] Files
  Options:
    --any-overlap-filter
      If not empty, set this FILTER if any variant in gnomad is found 
      overlaping the variant BUT we didn't find a correct match
      Default: <empty string>
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --bnd-distance
      Two BND variants are the same if their bounds are distant by less than 
      xxx bases. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 100
    --check-bnd-mate
      When comparing two BND, check that their mate (using the ALT allele) are 
      the same too
      Default: false
    --disable-annotations, -A
      do not transfer annotations from GNOMAD to VCF. Just a way to filter the 
      VCF using gnomad Frequency data
      Default: false
    --discordant_svtype
      If not empty, set this FILTER if SVTYPE are discordants
      Default: <empty string>
    --filter
      set this FILTER is the allele frequency found in the population is not 
      min-af<=x<=max-af. Discard variant if it is blank.
      Default: BAD_AF
    --force-svtype
      When comparing two SV variants, their INFO/SVTYPE should be the same. 
      Default is to just use coordinates to compare non-BND variants.
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -g, --gnomad
      Gnomad-SV VCF file. see 
      https://gnomad.broadinstitute.org/downloads#structural-variants 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --in-gnomad-filter
      If not empty, set this FILTER is variant was found in gnomad
      Default: <empty string>
    --max-af
      max allele frequency in watched population. A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
      Default: 1.0
    --min-af
      min allele frequency in watched population. A decimal number between 0.0 
      and 1.0. If the value ends with '%' it is interpretted as a percentage 
      eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' 
      => '0.01'.
      Default: 0.0
    -o, --out
      Output file. Optional . Default: stdout
    --population
      Watch gnomad population for AF
      Default: POPMAX_AF
    -p, --prefix
      INFO field prefix
      Default: GNOMAD_
    --sv-alleles-bases
      When comparing two non-BND SV variants, use their ALT alleles to adjust 
      the interval. It solves the problem of  
      'chr2:10556028:AATTATATAATTAAATTAATTATATAATT:A'  vs 
      'chr2:10556028:A:AATTATATAATTAAATTAATTATATAATT'. See 
      https://twitter.com/yokofakun/status/1169182491606999046 
      Default: false
    --sv-fraction
      Two SV have are the same if they share a fraction 'x' of their bases. 
      For very small SV the fraction can be quite small while for large SV the 
      fraction should be close to 1. The Syntax is the following : 
      (<MAX_SIZE_INCLUSIVE>:<FRACTION as double or percent>;)+ . For example 
      if the SV as a size of 99bp, the fraction used with be 0.6 for 
      '10:0.1;100:0.6;1000:0.9'. For the smallest size, a simple overlap is a 
      positive match.
      Default: 10:0.5;100:0.75;1000:0.8;10000:0.9
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad
 * sv



## Creation Date

20190814

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomadsv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Example:

```
java -jar dist/vcfgnomadsv.jar \
	-g src/test/resources/gnomad_v2_sv.sites.vcf.gz \
	./src/test/resources/manta.B00GWGD.vcf.gz
```


