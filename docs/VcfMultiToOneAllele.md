# VcfMultiToOneAllele

'one variant with N ALT alleles' to 'N variants with one ALT'


## Usage

```
Usage: vcfmulti2oneallele [options] Files
  Options:
    --addNoVariant
      Print Variants without ALT allele
      Default: false
    --disableHomVarAlt
      by default is a genotype is homvar for an external ALT ('2/2'), it will 
      be set to ./. (no call). Setting this option will replace the current 
      allele. 
      Default: true
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfmulti2oneallele
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAllele.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfMultiToOneAllele.java)


<details>
<summary>Git History</summary>

```
Fri Sep 8 12:42:11 2017 +0200 ; gnomad spring + add test ; https://github.com/lindenb/jvarkit/commit/03445831f08a7e61c34d0c6fab5c4c6b4d647c6c
Fri Sep 8 11:29:13 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/0c5035b59124ca18aba0405d0c616b565a32d10e
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Mon Jun 26 17:29:03 2017 +0200 ; burden ; https://github.com/lindenb/jvarkit/commit/a3b7abf21d07f0366e81816ebbb2cce26b2341e7
Fri Jun 23 16:37:51 2017 +0200 ; alt vs homref ; https://github.com/lindenb/jvarkit/commit/9ef5f8c8d0b33994515b0faac60e84b275ab34eb
Fri Jun 23 15:32:55 2017 +0200 ; updated vcf2multiallele ; https://github.com/lindenb/jvarkit/commit/9a69f3da7748ae43458e34f0cd3c0f052aa09b51
Fri Jun 23 15:26:55 2017 +0200 ; updated vcf2multiallele ; https://github.com/lindenb/jvarkit/commit/775e8ddcc38a3e283cf49d9287b06510d7634e31
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Thu Apr 27 17:22:22 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/0a27a246a537d2b48201596067652ea26bfc28d6
Fri Nov 25 12:30:34 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/31949a5be3c9948eb6d6fa72a96e8cbcbc66796d
Mon Jan 25 18:43:22 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/f49222d5a112d79fc4148f2a2a56e46a7ee5f517
Sat Jan 23 14:55:35 2016 +0100 ; factory builder ; https://github.com/lindenb/jvarkit/commit/d70912b7dbbca748cf4d45a0ba44a6bc70f804d7
Fri Nov 27 15:22:25 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/04a83d5b9f0e69fd2f7087e519b0de3e2b4f9863
Tue Jul 21 17:15:14 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc3230eabbfb7c2c9763528c63c1f42ae1281351
Mon Jul 6 16:14:07 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ee95fe6971b5655c61d7feb22e8fa877201a9ca6
Wed Jul 1 20:06:47 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/a98d085769cbd9a4300329cda346ab0697d24c61
Wed Jul 1 19:37:31 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/0c182acaac09c876387bbb4d0777fd6596284665
Tue Jun 30 17:45:37 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7c0d31b60217243a99bc4e1ea0045c9f885ba9bd
Mon Jun 29 16:03:33 2015 +0200 ; VcfMultiToOneAllele support samples in header ; https://github.com/lindenb/jvarkit/commit/70f05f90c0062763d380921d0c3a9cbcf33421a6
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Tue Jun 2 16:04:16 2015 +0200 ; 'one variant, N ALT' to 'N variants one ALT' #tweet ; https://github.com/lindenb/jvarkit/commit/6019125b77027d0d11ec86c6ff8de72413be7263
```

</details>

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


