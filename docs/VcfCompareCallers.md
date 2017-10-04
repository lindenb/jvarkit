# VcfCompareCallers

Compare two VCFs and print common/exclusive information for each sample/genotype


## Usage

```
Usage: vcfcomparecallers [options] Files
  Options:
    -B, --bed
      Limit to variants in that BED region
    --collapseGenotypeType
      collapse Genotype Type. Just show Same Genotype or Discordant, don't 
      print the type of genotype.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --jexl1
      An optional list of GATK-like JEXL expressions to filter the variants 
      from VCF File 1
      Default: []
    --jexl2
      An optional list of GATK-like JEXL expressions to filter the variants 
      from VCF File 2
      Default: []
    -c, --noCallIsHomRef
      No Call is HomRef (created when comparing merged vcf with GATK: there is 
      no homref, everything is nocall)
      Default: false
  * -o, --output
      Directory or zip file to save results to be plotted with gnuplot
    -p, --prefix
      Archive prefix (for option -d)
      Default: <empty string>
    -vcf1, --vcf1
      short descriptive name for VCF1
      Default: VCF1
    -vcf2, --vcf2
      short descriptive name for VCF2
      Default: VCF2
    --version
      print version and exit

```


## Keywords

 * vcf
 * compare
 * genotype


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
$ make vcfcomparecallers
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallers.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallers.java)


<details>
<summary>Git History</summary>

```
Thu Jul 13 20:16:36 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/85b6c9c196e9a065dfd47bee37fe50238af41660
Wed Jul 5 11:08:10 2017 +0200 ; vcffilterjdk ; https://github.com/lindenb/jvarkit/commit/b25cc45aa3a057f0dad46f0d83669bc88cc95e0c
Tue Jul 4 14:44:16 2017 +0200 ; rewritting vcfcomparecallers + jexl ; https://github.com/lindenb/jvarkit/commit/f8ec122b6f76218a66a7c8e7f7d5f4c203b56a9f
Tue Jul 4 12:30:21 2017 +0200 ; rewritting vcfcomparecallers ; https://github.com/lindenb/jvarkit/commit/b368bd58ec54e47911197870b82e7d78636aa3cf
Mon Jul 3 18:24:52 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/143424e4b26825a00eb8ac652f2b80ebe1ac79a8
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Wed Apr 26 17:26:23 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/ab6c7b760cd5376e08da24426cede7f84a6b3ae2
Fri Apr 21 18:16:07 2017 +0200 ; scan sv ; https://github.com/lindenb/jvarkit/commit/49b99018811ea6a624e3df556627ebdbf3f16eab
Wed Jan 18 18:26:00 2017 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/4f808db2ad9a8afff2f2e4bfc59cf4f11c29e0f9
Tue Feb 2 10:19:43 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/942c2c8b286da4d356535758ae16a8959c2ed58f
Wed Jan 13 15:25:58 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/db4a0f749e0c5b5a0ba067c7f4e89392ea6b62c3
Thu Nov 26 09:11:28 2015 +0100 ; continue central ; https://github.com/lindenb/jvarkit/commit/93af6fca2b99af07acf0216f783b270fd84dcaea
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Thu Mar 26 14:28:15 2015 +0100 ; xml output for VcfCompareCallers ; https://github.com/lindenb/jvarkit/commit/19dc404ebd0767a59897c510f900471a6ca2a47a
Wed Mar 25 16:04:36 2015 +0100 ; vcf-compare-callers: added more categories #tweet ; https://github.com/lindenb/jvarkit/commit/65f60dc554a14d7d3c72b5bc8cbfff416ac580b3
Tue Mar 24 17:42:28 2015 +0100 ; 1st commit for a tool comparing two vcfs (same samples but != callers)  #tweet ; https://github.com/lindenb/jvarkit/commit/7c64142e1df8f1ec072a1835386ad3b1f8fa237c
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparecallers** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Synopsis

```
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) stdin 
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) file2.vcf(.gz) 

```



both vcf **must** share the same sequence dictionary and must be sorted


#### History

* 20170704 : rewritten from scratch


### Example



```
$ java -jar dist/vcfcomparecallers.jar -o tmp Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz
$ (cd tmp && make)
```



