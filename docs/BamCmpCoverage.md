# BamCmpCoverage

Creates the figure of a comparative view of the depths sample vs sample. Memory consideration: the tool alloc an array of bits which size is: (MIN(maxdepth-mindepth,pixel_width_for_one_sample) * count_samples)^2


## Usage

```
Usage: bamcmpcoverage [options] Files
  Options:
    -b, --bed
      restrict to region
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -M, --maxDepth
      max depth
      Default: 1000
    -m, --minDepth
      min depth
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    -r, --region
      restrict to region
    --version
      print version and exit
    -w, --width
      image width
      Default: 1000

```


## Keywords

 * sam
 * bam
 * visualization
 * coverage


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
$ make bamcmpcoverage
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamCmpCoverage.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamCmpCoverage.java)


<details>
<summary>Git History</summary>

```
Wed Jun 21 17:31:49 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/556a9b4ed5d047a215e160c0a480ea241cea83d9
Thu Jun 15 15:30:26 2017 +0200 ; update vcfcalledwithanothermethod, vcfucsc ; https://github.com/lindenb/jvarkit/commit/0efbf47c1a7be8ee9b0a6e2e1dbfd82ae0f8508f
Mon May 15 10:41:51 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c13a658b2ed3bc5dd6ade57190e1dab05bf70612
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Fri May 5 15:06:21 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/4d2bbfed84609bdf14eb1b14a35ab24eb8ad5b26
Fri Jun 17 13:56:39 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/865252a44fc018f46b4280788cec65a1383dcc18
Tue Nov 3 22:42:18 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/4e4a9319be20626f0ea01dc2316c6420ba8e7dac
Wed Sep 23 18:01:13 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/6156d1359a63a80c24f5b7694dc70431f6816289
Tue Jun 9 12:17:32 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3601851f8d35e25d0130b1cb765c936e53292750
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Wed Dec 3 12:56:00 2014 +0100 ; first hilbert, vcf detect fileformat, label height dans bamcmpdepth ; https://github.com/lindenb/jvarkit/commit/eb3aaf4591d7f520438521edf07751ba6968731c
Tue Dec 2 14:35:01 2014 +0100 ; side effec in bamcmpcov ; https://github.com/lindenb/jvarkit/commit/8fb5e7542d432f8d1d15a49a03fea4c3c0ccc501
Tue Dec 2 12:41:22 2014 +0100 ; improved bamcmpcov ; https://github.com/lindenb/jvarkit/commit/a54af7812f0c12492ea4824ac25ea6b7bf92f33d
Fri Nov 28 17:07:19 2014 +0100 ; binning bam cov ; https://github.com/lindenb/jvarkit/commit/3c2fd832f2f2402638b50d0f02759099f00ca048
Thu Nov 27 16:41:49 2014 +0100 ; bam cmp coverage ; https://github.com/lindenb/jvarkit/commit/a820811273daf971c0baceaf63778934a43da070
Thu Nov 27 13:11:06 2014 +0100 ; bam compare coverage ; https://github.com/lindenb/jvarkit/commit/0be60cca2b40fa2bb2713e759271573936911aba
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamcmpcoverage** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)


```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```






### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)

### Example

```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```



