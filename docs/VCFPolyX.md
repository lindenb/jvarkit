# VCFPolyX

Number of repeated REF bases around POS.


## Usage

```
Usage: vcfpolyx [options] Files
  Options:
    -n, --filter
      if number of repeated bases is greater or equal to 'n' set a FILTER = 
      (tag) 
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -t, --tag
      Tag used in INFO and FILTER columns.
      Default: POLYX
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
 * repeat


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
$ make vcfpolyx
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyX.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VCFPolyX.java)

Git History for this file:
```
Mon Aug 7 09:53:19 2017 +0200 ; fixed unicode problems after https://github.com/lindenb/jvarkit/issues/82 ; https://github.com/lindenb/jvarkit/commit/68254c69b027a9ce81d8b211447f1c0bf02dc626
Thu Jul 13 18:38:11 2017 +0200 ; bioalcidaejdk / vcfstats ; https://github.com/lindenb/jvarkit/commit/b403134b5c8e9961489ac8b41d477947a83ff2c4
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Sat Apr 29 18:45:47 2017 +0200 ; partition ; https://github.com/lindenb/jvarkit/commit/7d72633d50ee333fcad0eca8aaa8eec1a475cc4d
Fri Apr 14 17:09:04 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/881f52d5d3775325240114702b7f07148b626f4c
Tue Mar 22 17:19:22 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/97e0e23bddd49049c71d56d495d090c0af636670
Fri Dec 4 12:12:44 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/bfa0325927a7be363feecc143ef1d8a3de71f483
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Thu Apr 30 16:54:43 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/eb6a52869eb3b9b8bf048a079e1d7a96bccab2bf
Wed Apr 29 09:42:54 2015 +0200 ; moved tools vcfpolyx & vcfbed to std argc/argv #tweet ; https://github.com/lindenb/jvarkit/commit/84c8fcb3aa33e74c0a92c6986629a761330c9afe
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Fri Oct 25 17:42:45 2013 +0200 ; close Reference Fasta (picard.100) ; https://github.com/lindenb/jvarkit/commit/9c4a6831016175308ec9a80539e2093c32e78af9
Fri Oct 11 15:39:02 2013 +0200 ; picard v.100: deletion of VcfIterator :-( ; https://github.com/lindenb/jvarkit/commit/e88fab449b04aed40c2ff7f9d0cf8c8b6ab14a31
Mon Sep 30 17:03:13 2013 +0200 ; VCFPoly-X added ; https://github.com/lindenb/jvarkit/commit/9c91e722648f8180b964d52579c326fe583b010c
```

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpolyx** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java  -jar dist/vcfpolyx.jar -R reference.fa input.vcf
(...)
2   1133956 .   A   G   2468.84 .   POLYX=23
2   1133956 .   A   AG  3604.25 .   POLYX=23
2   2981671 .   T   G   47.18   .   POLYX=24
(...)
```


