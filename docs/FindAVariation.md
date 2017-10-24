# FindAVariation

Finds a specific mutation in a list of VCF files


## Usage

```
Usage: findavariation [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -homref, --homref
      Hide HOM_REF genotypes
      Default: false
    -indexed, --indexed
      [20171020] Search only in indexed vcf
      Default: false
    -nocall, --nocall
      Hide NO_CALL genotypes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -f, --posfile
      Add this file containing chrom:position
      Default: []
    -p, --position
      A list of 'chrom/position'
      Default: []
    -snp, --snp
      Search only variant have the very same position (ignore overlapping 
      variants) 
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * variation
 * search


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
$ make findavariation
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java)


<details>
<summary>Git History</summary>

```
Fri Oct 20 09:53:19 2017 +0200 ; skat continue ; https://github.com/lindenb/jvarkit/commit/76b0e511e054e438c38d2157bbc0e148480288bb
Mon Sep 11 14:48:00 2017 +0200 ; adding tests, add test files for gnomad ; https://github.com/lindenb/jvarkit/commit/bc90c3c76e38e677a2fe824ce29bd7705dde3bd0
Thu Jul 13 20:16:36 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/85b6c9c196e9a065dfd47bee37fe50238af41660
Mon May 15 10:41:51 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c13a658b2ed3bc5dd6ade57190e1dab05bf70612
Wed Apr 5 13:49:50 2017 +0200 ; cont, fix bug in findallcovatpos ; https://github.com/lindenb/jvarkit/commit/7db18c7fe90fd5bf64d3ff3a4505607a1974ce6b
Fri Mar 31 17:08:11 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/f78937d19c4b038e69a32fbcfa2aeab8fd8417c6
Thu Mar 30 17:38:36 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/bba625df69e00a0aa54de192cdce6fda110a65b4
Tue Mar 22 17:19:22 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/97e0e23bddd49049c71d56d495d090c0af636670
Fri Nov 27 15:22:25 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/04a83d5b9f0e69fd2f7087e519b0de3e2b4f9863
Tue Nov 24 16:06:19 2015 +0100 ; fix https://github.com/lindenb/jvarkit/issues/36 ; https://github.com/lindenb/jvarkit/commit/eac04e587d9e0f784dd1a00c2d1245891a537568
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Fri Nov 28 12:44:44 2014 +0100 ; find all coverages ; https://github.com/lindenb/jvarkit/commit/a8c96e489787bf94d752e6bbd7c091175617459b
Mon Jun 23 17:34:39 2014 +0200 ; jnlp ; https://github.com/lindenb/jvarkit/commit/fbd79898e088a286a8ccb474ca0eed1a7d64d876
Mon Jun 23 12:34:44 2014 +0200 ; find-a-variation + using abstractcodec instead of vcfcodec ; https://github.com/lindenb/jvarkit/commit/da621ba8326d56da8f6907c845c539e4ea785284
```

</details>

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
   java -jar dist/findamutation.jar -p "chr1:1234" 


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

 

