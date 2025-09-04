# SVCasesControls

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find SV present in cases but not in controls.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar svcasescontrols  [options] Files

Usage: svcasescontrols [options] Files
  Options:
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
    --complex
      By default this tool select the SV that are found in the cases but not 
      in the controls. When using this flag, all variants for cases are 
      extracted and a count of CASE/CONTROL having the SV is added in the INFO 
      column. 
      Default: false
    -c, --contig, --chrom
      limit to this contig
    --force-svtype
      When comparing two SV variants, their INFO/SVTYPE should be the same. 
      Default is to just use coordinates to compare non-BND variants.
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --no-bnd
      discar BND
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
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

 * sv
 * manta
 * vcf



## Creation Date

20240513

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SVCasesControls.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SVCasesControls.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **svcasescontrols** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
# Motivation
 
Find SV present in cases and controls.


 
# Input
 
 input is a tab delimited file with a header:
 
  - vcf: the path to the vcf
  - sample : the sample name associated to the vcf. If empty the name will be searched in the CHROM header
  - status : case OR control
 
 
# Example
 
 ```

 $ java -jar jvarkit.jar svcasescontrols input.tsv > output.vcf
 ```


