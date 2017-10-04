# VCFCompareGT

 compare two or more genotype-callers for the same individuals. Produce a VCF with FORMAT fields indicating if a genotype is new or modified.


## Usage

```
Usage: vcfcomparegt [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -m
      only print modified samples
      Default: false

```


## Keywords

 * vcf
 * compare


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
$ make vcfcomparegt
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGT.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGT.java)


<details>
<summary>Git History</summary>

```
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Fri Apr 21 18:16:07 2017 +0200 ; scan sv ; https://github.com/lindenb/jvarkit/commit/49b99018811ea6a624e3df556627ebdbf3f16eab
Fri Jan 22 23:49:23 2016 +0100 ; vcfiterator is now an interface ; https://github.com/lindenb/jvarkit/commit/9f9b9314c4b31b21044c5911a7e79e1b3fb0af7a
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Tue Feb 24 16:43:03 2015 +0100 ; vcfin : code rewrittern. picky with ALT alleles. #tweet ; https://github.com/lindenb/jvarkit/commit/65ef7741539e89c7a1a1f9cca28c13d531902c96
Thu Sep 11 09:36:01 2014 +0200 ; problem with java dataInputSTream: writeUTF requires line.length < SHORt_MAX ; https://github.com/lindenb/jvarkit/commit/19eac4ee36909a730903546b50461de3c19a5c1f
Mon May 12 15:27:08 2014 +0200 ; moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/fd30a81154a16835b5bab3d8e1ef90c9fee6bdcb
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Wed Jan 29 10:50:00 2014 +0100 ; world map ; https://github.com/lindenb/jvarkit/commit/3fb0f8ad813d25ee0871e7e24c42693e1036438f
Wed Nov 6 06:25:34 2013 +0100 ; fix vcfcmp ; https://github.com/lindenb/jvarkit/commit/10808a65b7d25b34072d636912bf19f9edd0556f
Tue Nov 5 13:17:33 2013 +0100 ; vcf compare genotype ; https://github.com/lindenb/jvarkit/commit/89111e42853993ae28e35ea966e2c1d7306d1e57
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparegt** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash

$ java -jar dist/vcfcomparegt.jar -m  Sample.samtools.vcf.gz Sample.gatk.vcf.gz

##fileformat=VCFv4.1
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GCH,Number=1,Type=Integer,Description="Changed Genotype">
##FORMAT=<ID=GNW,Number=1,Type=Integer,Description="Genotype Created/Deleted">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Qual">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=GDF,Number=.,Type=String,Description="Samples with Genotype Difference">
##VCFCompareGT_1=File: Sample.samtools.vcf.gz
##VCFCompareGT_2=File: Sample.gatk.vcf.gz
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_2	Sample_1
X	1860854	rs5781	A	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1866893	rs2824	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1878904	.	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	0/1:20:0:1:71	./.
X	1895117	.	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/0:2:0:1:27
X	1895755	.	C	AG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:4:0:1:17
X	1900009	rs6181	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:13:0:1:30	./.
X	1905130	.	AG	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:3:0:1:16
X	1905160	.	A	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1905165	.	C	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
X	1913889	.	C	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1948846	rs6	T	TG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:239:0:1:99	./.
X	1955199	.	C	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
(...)
```


