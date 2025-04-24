# VcfUkbiobank

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

annotates an VCF with the https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!)


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfukbb  [options] Files

Usage: vcfukbb [options] Files
  Options:
    --api
      API base url
      Default: https://afb.ukbiobank.ac.uk/api
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --debug
      debug
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
    --retry
      max-retry
      Default: 10
    --seconds
      wait 'x' seconds between each call to the API
      Default: 10
    --skip-filtered
      Skip FILTERed variants
      Default: false
    --version
      print version and exit

```


## Keywords

 * ukbiobank
 * vcf



## Creation Date

20250424

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfukbb/VcfUkbiobank.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfukbb/VcfUkbiobank.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfukbb** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





# WARNING

I cannot explain why this software doesn't run with jdk17 (I always get a http error 403) but runs with jdk23...

The server might not like too many requests. Use a your own risk.

## Example:

```
$ cat jeter.vcf | grep -v "#"
chr22	10510356	.	T	A	.	.	.
chr22	10526438	.	A	G	.	.	.
```

```
java -jar jvarkit.jar vcfukbb jeter.vcf | grep -v "#"
chr22	10510356	.	T	A	.	.	UKBB_AFR_AC=288;UKBB_AFR_AF=0.319;UKBB_ALL_AC=12606;UKBB_ALL_AF=0.915;UKBB_ASJ_AC=71;UKBB_ASJ_AF=0.934;UKBB_EAS_AC=779;UKBB_EAS_AF=0.986;UKBB_NFE_AC=9989;UKBB_NFE_AF=0.954;UKBB_OTH_AC=448;UKBB_OTH_AF=0.907;UKBB_SAS_AC=1031;UKBB_SAS_AF=0.982
chr22	10526438	.	A	G	.	.	UKBB_AFR_AC=13052;UKBB_AFR_AF=0.863;UKBB_ALL_AC=536812;UKBB_ALL_AF=0.986;UKBB_ASJ_AC=3261;UKBB_ASJ_AF=0.964;UKBB_EAS_AC=3119;UKBB_EAS_AF=0.960;UKBB_NFE_AC=495657;UKBB_NFE_AF=0.990;UKBB_OTH_AC=9134;UKBB_OTH_AF=0.964;UKBB_SAS_AC=12589;UKBB_SAS_AF=0.991
```

##


