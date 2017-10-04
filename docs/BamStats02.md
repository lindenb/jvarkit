# BamStats02

Statistics about the flags and reads in a BAM


## Usage

```
Usage: bamstats02 [options] Files
  Options:
    -B, --bed
      Optional Bed File
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```

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
$ make bamstats02
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats02.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats02.java)


<details>
<summary>Git History</summary>

```
Fri Jun 2 16:31:30 2017 +0200 ; circos / lumpy ; https://github.com/lindenb/jvarkit/commit/7bddffca3899196e568fb5e1a479300c0038f74f
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Thu Apr 20 19:19:21 2017 +0200 ; fix ; https://github.com/lindenb/jvarkit/commit/793db1bdad0bcff729c6c5699d7623d44110eb02
Thu Apr 20 17:17:22 2017 +0200 ; continue transition jcommander ; https://github.com/lindenb/jvarkit/commit/fcf5def101925bea9ddd001d8260cf65aa52d6a0
Wed Mar 23 18:28:19 2016 +0100 ; pack ; https://github.com/lindenb/jvarkit/commit/f4e198c21581c9bd26d1c844122db75771e54cbd
Mon Apr 13 10:32:12 2015 +0200 ; added option -n for https://github.com/lindenb/jvarkit/issues/26#issuecomment-92193999 ; https://github.com/lindenb/jvarkit/commit/19072426981ab0f13e61755012a4a35da591cc95
Thu Apr 2 16:18:39 2015 +0200 ; read-name parser, enhanced bamstats02 with lane, flowcell, run etc... #tweet ; https://github.com/lindenb/jvarkit/commit/10834a090bb63fd6d0870eecd2cae8bd0211062a
Wed Apr 1 12:04:06 2015 +0200 ; bamstats02: statistics about the reads in one or more BAM #tweet ; https://github.com/lindenb/jvarkit/commit/b79a5374acce8261f249739c8fdd91e19ff4dac4
Tue Mar 31 22:04:56 2015 +0200 ; GUI exploring samflags/sample/chrom/mapq #tweet ; https://github.com/lindenb/jvarkit/commit/919a2eff175fddf74ea981d02de85813636596f2
Tue Mar 31 19:05:54 2015 +0200 ; bamstats02 ; https://github.com/lindenb/jvarkit/commit/3f01c89c98c2df5be453f2dd5f16d5eb4cced6a4
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamstats02** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$  find dir -name "*final.bam" | xargs  java -jar dist/bamstats02.jar -B capture.bed  > output.tsv
$  verticalize output.tsv

>>> 2
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      2
$4      mapq    60
$5      inTarget        1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      1
$8      READ_UNMAPPED   0
$9      MATE_UNMAPPED   0
$10     READ_REVERSE_STRAND     1
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   0
$13     SECOND_IN_PAIR  1
$14     NOT_PRIMARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   463982
<<< 2

>>> 3


>>> 3
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      .
$4      mapq    0
$5      inTarget        -1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      0
$8      READ_UNMAPPED   1
$9      MATE_UNMAPPED   1
$10     READ_REVERSE_STRAND     0
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   1
$13     SECOND_IN_PAIR  0
$14     NOT_PRIMARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   458630
<<< 3

>>> 4
$1      #filename       dir/Sample0258.final.bam
$2      sampleName      Sample0258
$3      chromosome      .
$4      mapq    0
$5      inTarget        -1
$6      READ_PAIRED     1
$7      READ_MAPPED_IN_PROPER_PAIR      0
$8      READ_UNMAPPED   1
$9      MATE_UNMAPPED   1
$10     READ_REVERSE_STRAND     0
$11     MATE_REVERSE_STRAND     0
$12     FIRST_IN_PAIR   0
$13     SECOND_IN_PAIR  1
$14     NOT_PRIMARY_ALIGNMENT   0
$15     READ_FAILS_VENDOR_QUALITY_CHECK 0
$16     READ_IS_DUPLICATE       0
$17     SUPPLEMENTARY_ALIGNMENT 0
$18     count   458630
<<< 4
```
```





### See also

BamStats02View



