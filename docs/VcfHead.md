# VcfHead

print the first variants of a vcf


## Usage

```
Usage: vcfhead [options] Files
  Options:
    -c, --bycontig
      number of variants
      Default: false
    -n, --count
      number of variants
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfhead
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfHead.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfHead.java)

Git History for this file:
```
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Fri Apr 14 15:27:32 2017 +0200 ; annotation proc ; https://github.com/lindenb/jvarkit/commit/72b9383a8472e5a91120bab84d15b8acad4db8d4
Mon Mar 27 17:55:22 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3b7033dcb365493f62d7b78ca4410b6bf3cd716d
Mon Mar 27 09:23:26 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/97ebb85e7afc877ed11317f2078e48104983a50c
Mon Mar 27 08:21:34 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/f65f2fab10d0d3e895f85a5d5fa433699d64d48f
Fri Apr 29 17:26:25 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/231856146f1035e21b4adfbc4e87a01b60d0d39e
Wed Feb 10 10:00:55 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/0a518250fb9b2e9894b9329d01afe1dfe0a2f6de
Tue Nov 3 22:42:18 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/4e4a9319be20626f0ea01dc2316c6420ba8e7dac
Wed Sep 23 15:36:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/03cc568ed51cf43d7af6f0d913a0f8c52ddfe5d7
Wed Sep 23 13:19:34 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/58308bb7b5859c5d3360a33e22f4d62607997a7c
Thu Mar 12 16:57:07 2015 +0100 ; tool to compare VCF with one sample called with multiple methods #tweet ; https://github.com/lindenb/jvarkit/commit/351c259dc9f1d8bebab19b3dc57fc6a610257542
Fri Feb 20 17:50:08 2015 +0100 ; continue integration in knime ; https://github.com/lindenb/jvarkit/commit/5693e07342a30f21c807a4e3a655e3446019458f
Fri Feb 20 13:25:12 2015 +0100 ; moving vcfhead, vcftail to knime. samsequenceprogress compact watch form ; https://github.com/lindenb/jvarkit/commit/34c137e5f12c05b7eb9d7ce1b546cde6c8890cc3
Fri Nov 7 17:11:35 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/17b77998c42387988386f5a696bba464d130cf86
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Tue Dec 10 22:12:54 2013 +0100 ; vcf head tail ; https://github.com/lindenb/jvarkit/commit/25adbeaf3de1c64fd793c5ddf5236322adb430d1
```

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfhead** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 

```bash
$ curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
 java -jar dist/vcfhead.jar -n 2 | grep -v "##"

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  1094PC0005  1094PC0009  1094PC0012  1094PC0013
chr1    30860   .   G   C   33.46   .   AC=2;AF=0.053;AN=38;BaseQRankSum=2.327;DP=49;Dels=0.00;EFF=DOWNSTREAM(MODIFIER||||85|FAM138A|protein_coding|COD
ING|ENST00000417324|),DOWNSTREAM(MODIFIER|||||FAM138A|processed_transcript|CODING|ENST00000461467|),DOWNSTREAM(MODIFIER|||||MIR1302-10|miRNA|NON_CODING
|ENST00000408384|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000469289|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST000004
73358|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING
|ENST00000430492|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene
|NON_CODING|ENST00000488147|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000538476|);FS=3.128;HRun=0;HaplotypeScore=0.6718;In
breedingCoeff=0.1005;MQ=36.55;MQ0=0;MQRankSum=0.217;QD=16.73;ReadPosRankSum=2.017 GT:AD:DP:GQ:PL  0/0:7,0:7:15.04:0,15,177    0/0:2,0:2:3.01:0,3,39   0
/0:6,0:6:12.02:0,12,143    0/0:4,0:4:9.03:0,9,119
chr1    69270   .   A   G   2694.18 .   AC=40;AF=1.000;AN=40;DP=83;Dels=0.00;EFF=SYNONYMOUS_CODING(LOW|SILENT|tcA/tcG|S60|305|OR4F5|protein_coding|CODI
NG|ENST00000335137|exon_1_69091_70008);FS=0.000;HRun=0;HaplotypeScore=0.0000;InbreedingCoeff=-0.0598;MQ=31.06;MQ0=0;QD=32.86 GT:AD:DP:GQ:PL  ./. ./. 1/
1:0,3:3:9.03:106,9,0  1/1:0,6:6:18.05:203,18,0
```
 

