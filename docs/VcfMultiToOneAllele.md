# VcfMultiToOneAllele

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

'one variant with N ALT alleles' to 'N variants with one ALT'


## Usage

```
Usage: vcfmulti2oneallele [options] Files
  Options:
    --addNoVariant
      Print Variants without ALT allele
      Default: false
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    --disableHomVarAlt
      by default is a genotype is homvar for an external ALT ('2/2'), it will 
      be set to ./. (no call). Setting this option will replace the current 
      allele. 
      Default: true
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -highest, --highest
      [20170723]: Use  Allele With Highest Allele Count, discard/replace the 
      other 
      Default: false
    --ignoreMissingInfoDecl
      Ignore error when a variant INFO is missing a definition in the VCF 
      header. 
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    --replaceWith
      When replacing an alternative allele, replace it with REF or current ALT 
      allele. 
      Default: REF
      Possible Values: [REF, ALT]
    -r, --rmAtt
      [20161110]: after merging with GATK CombineVariants there can have 
      problemes with INFO/type='A' present in vcf1 but not in vcf2, and 
      multiallelelic variants. This option delete the attributes having such 
      problems. 
      Default: false
    -p, --samples
      print sample genotypes.
      Default: false
    --skipSpanningDeletions
      Skip Alt Spanning deletion alleles *
      Default: false
    -tag, --tag
      Info field name that will be added to recall the original alleles.
      Default: VCF_MULTIALLELIC_SRC
    --vc-attribute-recalc-ignore-filtered
      When recalculating variant attributes like DP AF, AC, AN, ignore 
      FILTERed **Genotypes**
      Default: false
    --vc-attribute-recalc-ignore-missing
      Ignore missing VCF headers (DP, AF, AC, AN). Default behavior: adding 
      VCF header if they're missing
      Default: false
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfmulti2oneallele
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAllele.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAllele.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAlleleTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAlleleTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmulti2oneallele** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## SNPEFF/VEP Annotations

this tool will try to split the VEP or the VCF annotation for each allele.


## Example

Exac contains multi-ALT  variants:

```bash
$ gunzip -c ExAC.r0.3.sites.vep.vcf.gz | grep rs3828049

1	889238	rs3828049	G	A,C	8422863.10	PASS	AC=6926,3;AC_AFR=220,0;AC_AMR=485,1;AC_Adj=6890,3;AC_EAS=746,0;AC_FIN=259,0;AC_Het=6442,3,0;AC_Hom=224,0;AC_NFE=3856,0;AC_OTH=41,0;AC_SAS=1283,2;AF=0.057,2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;(...)

```

processed with this tools:
```
$ java -jar dist/vcfmulti2oneallele.jar  ExAC.r0.3.sites.vep.vcf.gz   | grep rs3828049

1	889238	rs3828049	G	A	8422863.10	PASS	AC=6926;AC_AFR=220;AC_AMR=485;AC_Adj=6890;AC_EAS=746;AC_FIN=259;AC_Het=6442;AC_Hom=224;AC_NFE=3856;AC_OTH=41;AC_SAS=1283;AF=0.057;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;BaseQRankSum=-2.170e-01;VCF_MULTIALLELIC_SRC=A|C;(...)
1	889238	rs3828049	G	C	8422863.10	PASS	AC=3;AC_AFR=0;AC_AMR=1;AC_Adj=3;AC_EAS=0;AC_FIN=0;AC_Het=3;AC_Hom=0;AC_NFE=0;AC_OTH=0;AC_SAS=2;AF=2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;VCF_MULTIALLELIC_SRC=A|C;(....)
```

## History

* 20170606 added support for VCFHeaderLineCount.R

