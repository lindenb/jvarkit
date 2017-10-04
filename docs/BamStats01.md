# BamStats01

Statistics about the reads in a BAM.


## Usage

```
Usage: samstats01 [options] Files
  Options:
    -B, --bed
      capture bed file. Optional
    --groupby
      Group Reads by
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -q, --qual
      min mapping quality
      Default: 30.0
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
$ make samstats01
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java)


<details>
<summary>Git History</summary>

```
Fri Jun 30 17:20:28 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/5fe0a8568d706e8cd5898ff66c0ebd8b1f8447a5
Fri Jun 2 16:31:30 2017 +0200 ; circos / lumpy ; https://github.com/lindenb/jvarkit/commit/7bddffca3899196e568fb5e1a479300c0038f74f
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Sat Apr 29 18:45:47 2017 +0200 ; partition ; https://github.com/lindenb/jvarkit/commit/7d72633d50ee333fcad0eca8aaa8eec1a475cc4d
Wed Apr 19 17:58:48 2017 +0200 ; rm-xml ; https://github.com/lindenb/jvarkit/commit/95f05cfd4e04f5013c22274c49db7bcc4cbbb1c8
Wed Jun 8 12:51:03 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3a139dad3aa0c899b4a84c9a0d2908d47ecccd58
Fri Apr 15 17:09:59 2016 +0200 ; inject pedigree ; https://github.com/lindenb/jvarkit/commit/f9a18a1ce155a78b2e430d8d7860d0cab8f33722
Fri Jan 23 09:07:46 2015 +0100 ; MAPQ=0 in samstats01, updated SequenceOntology , updated vep ; https://github.com/lindenb/jvarkit/commit/388c0adae3f111c3d77819a19c756755e40bf5f3
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Apr 2 17:50:50 2014 +0200 ; cont rnaseq ; https://github.com/lindenb/jvarkit/commit/7b3f7e13a112b09018284931678ac78dd32cefcc
Wed Dec 11 15:01:44 2013 +0100 ; improv bam-gui. New sam flag from bamstats01 ; https://github.com/lindenb/jvarkit/commit/314bf88924a4003e6d6189ad3280d8b4df485aa1
Tue Dec 10 12:30:02 2013 +0100 ; cont. ; https://github.com/lindenb/jvarkit/commit/88a74712aabfc418d7ed390fff6e965e7491b210
Mon Dec 9 11:37:46 2013 +0100 ; vcf2xml , bamstats01 X/Y ; https://github.com/lindenb/jvarkit/commit/2c13f6f369faf3d076ccc9420b5284cd990c6892
Fri Dec 6 15:55:48 2013 +0100 ; bamstats 01 : count X and Y ; https://github.com/lindenb/jvarkit/commit/6c29728c3b83cde2247659e1c39aa355971f1f6d
Wed Dec 4 17:08:33 2013 +0100 ; cgi ; https://github.com/lindenb/jvarkit/commit/0d5df7769ee9a585ffa51dce809070f3f222f02a
Thu Nov 7 13:54:08 2013 +0100 ; continue ; https://github.com/lindenb/jvarkit/commit/2a7db844ef92646208fb98090906fdb21163613d
Mon Nov 4 16:23:17 2013 +0100 ; moving to std cli program ; https://github.com/lindenb/jvarkit/commit/9ed1dc3f053d57d379862c9a14648f96a967ada7
Mon Nov 4 13:49:15 2013 +0100 ; chaned command line handling + getopt ; https://github.com/lindenb/jvarkit/commit/939d2ccf1a9a4be2d2116586b925062c65d81195
Fri Oct 11 15:39:02 2013 +0200 ; picard v.100: deletion of VcfIterator :-( ; https://github.com/lindenb/jvarkit/commit/e88fab449b04aed40c2ff7f9d0cf8c8b6ab14a31
Mon Sep 16 09:06:18 2013 +0200 ; bamstats01 ; https://github.com/lindenb/jvarkit/commit/f38435a64d5a4456e46c0e37ff1fbe115bfcfa3e
Mon Sep 16 00:12:39 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ca35f65cac4a74a3873a14ba036a6856b1fd02fc
Fri Jul 26 17:22:50 2013 +0200 ; vendredi 17H00 ... dodo ; https://github.com/lindenb/jvarkit/commit/347e4d3b00d3642542d9cca6723bd8157c395da6
Fri Jul 5 08:12:55 2013 +0200 ; misc code added ; https://github.com/lindenb/jvarkit/commit/73fae24c634a5275693493a3899449c78a2474e1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samstats01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### History



* Dec 2013 Added PROPER_PAIR_HMQ for @SolenaLS
* Dec 2013 Added X and Y for @SolenaLS



### Output


See also: http://picard.sourceforge.net/explain-flags.html


#### Counts


* TOTAL : total number of reads (not PAIRS of reads)
* PAIRED: total number of reads paired (should be equals to ALL for a paired-end assay)
* UNMAPPED : count unmapped reads 
* MAPPED  : count mapped reads
* PROPER_PAIR  : count reads in proper pair (forward+reverse, same chromosome, distance is OK)
* PROPER_PAIR_HMQ  : proper pairs with mapping quality >= user qual
* PLUS : reads on plus strand
* MINUS : reads on minus strand
* PRIMARY_ALIGNMENT : alignment flagged as primary alignment (not alternative position)
* FAIL_MAPPING_QUALITY : MAQ < user qual
* DUPLICATE : the flag 'duplicate' was set
* FAIL_VENDOR_QUALITY : the flag "read fails platform/vendor quality checks" was set
* OK_FOR_PE_CALLING : reads ok for Paired-end mapping ( properly paired, not dup, not fails_vendor_qual,  not fails_mapping_qual, primary align )
* X and Y : number of reads mapping the chromosomes X/chrX and Y/chrY


#### Categories


* ALL: all reads
* IN_TARGET: reads overlapping user's BED (if provided)
* OFF_TARGET: reads with no overlap with user's BED (if provided)



### Example


```
$  java -jar dist/bamstats01.jar \
		IN=my.bam \
		BED=capture.bed

(...)
#Filename	Sample	ALL_TOTAL	ALL_PAIRED	ALL_UNMAPPED	ALL_MAPPED	ALL_PROPER_PAIR	ALL_PLUS_STRAND	ALL_MINUS_STRAND	ALL_PRIMARY_ALIGNMENT	ALL_FAIL_MAPPING_QUALITY	ALL_DUPLICATE	ALL_FAIL_VENDOR_QUALITY	IN_TARGET_TOTAL	IN_TARGET_PAIRED	IN_TARGET_UNMAPPED	IN_TARGET_MAPPED	IN_TARGET_PROPER_PAIR	IN_TARGET_PLUS_STRAND	IN_TARGET_MINUS_STRAND	IN_TARGET_PRIMARY_ALIGNMENT	IN_TARGET_FAIL_MAPPING_QUALITY	IN_TARGET_DUPLICATE	IN_TARGET_FAIL_VENDOR_QUALITY	OFF_TARGET_TOTAL	OFF_TARGET_PAIRED	OFF_TARGET_UNMAPPED	OFF_TARGET_MAPPED	OFF_TARGET_PROPER_PAIR	OFF_TARGET_PLUS_STRAND	OFF_TARGET_MINUS_STRAND	OFF_TARGET_PRIMARY_ALIGNMENT	OFF_TARGET_FAIL_MAPPING_QUALITY	OFF_TARGET_DUPLICATE	OFF_TARGET_FAIL_VENDOR_QUALITY
my.bam	Sample	1617984	1617984	3966	1614018	1407862	806964	807054	1614018	56980	0	0	1293922	1293922	0	1293922	1133808	644741	649181	1293922	14087	0	0	320096	320096	0	320096	274054	162223	157873	320096	42893	0	0
(...)

```



