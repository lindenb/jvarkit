# VcfToBam

vcf to bam


## Usage

```
Usage: vcf2bam [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --fragmentsize
      fragment size
      Default: 600
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --readsize
      read size
      Default: 100
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --version
      print version and exit

```


## Keywords

 * ref
 * vcf
 * bam


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
$ make vcf2bam
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToBam.java)


<details>
<summary>Git History</summary>

```
Mon Sep 4 17:34:39 2017 +0200 ; fix https://github.com/lindenb/jvarkit/issues/86#issuecomment-326986654 ; https://github.com/lindenb/jvarkit/commit/7e9a296f733cfa76364b92b34707ec33d8e26f64
Sun Sep 3 00:12:21 2017 +0200 ; fix https://github.com/lindenb/jvarkit/issues/86 ; https://github.com/lindenb/jvarkit/commit/28ae7e722db261d7d337e066f52bfb9d88e53733
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Tue Apr 18 18:26:58 2017 +0200 ; which changes ?? ; https://github.com/lindenb/jvarkit/commit/2d7cf86faca95815601e4bdd516a757c960749a3
Fri Oct 2 18:46:06 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/470e305ccf3036229546d3f3232d5cc8b230fc27
Thu Jun 25 16:18:29 2015 +0200 ; extends REF sequence with clipped reads #tweet ; https://github.com/lindenb/jvarkit/commit/e3e4b7c31e357848b2e156affaaead86a8b5cefe
Fri Jun 12 21:12:27 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/56febbf1d9207f523d3ce342ca6c7b7ecf681fcc
Fri Jun 12 18:32:07 2015 +0200 ; starting vcf to bam ; https://github.com/lindenb/jvarkit/commit/dfa534e03f973083f41247bdae20637f6232a358
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2bam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

```bash
$ java -jar dist-1.139/vcf2bam.jar -R ref.fa  samtools.vcf.gz 2> /dev/null | grep -v "100="



@HD	VN:1.5	SO:unsorted
@SQ	SN:seg1	LN:5101
@SQ	SN:seg2	LN:4000
@RG	ID:NA12878	SM:NA12878	LB:illumina	DS:NA12878
@RG	ID:NA12891	SM:NA12891	LB:illumina	DS:NA12891
@RG	ID:NA12892	SM:NA12892	LB:illumina	DS:NA12892
@CO	Generated with -R /home/lindenb/src/ngsxml/test/ref/ref.fa /home/lindenb/src/ngsxml/OUT/Projects/Proj1/VCF/samtools/Proj1.samtools.vcf.gz
0000003609	83	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003610	147	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003611	147	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003612	83	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003615	83	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003616	147	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003617	147	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003627	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003629	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003630	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003633	147	seg1	1453	60	95=1X4=	=	953	-500	TTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003635	83	seg1	1453	60	95=1X4=	=	953	-500	TTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003640	147	seg1	1454	60	94=1X5=	=	954	-500	TAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTCA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1

```



