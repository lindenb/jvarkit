# VcfUcsc

annotate an VCF with mysql UCSC data


## Usage

```
Usage: vcfucsc [options] Files
  Options:
    -D, --database
      database name
      Default: hg19
    -e, --expression
      expression string.
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -T, -t, --table
      table name
    --version
      print version and exit

```


## Keywords

 * ucsc
 * mysql
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
$ make vcfucsc
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java)


<details>
<summary>Git History</summary>

```
Fri Jun 16 12:34:13 2017 +0200 ; xsltstream ; https://github.com/lindenb/jvarkit/commit/e94b85ebca08b23419359bf15e134a6f63823582
Fri Jun 16 09:26:28 2017 +0200 ; old ; https://github.com/lindenb/jvarkit/commit/a26808647296127b65d358c382c6aa6acf3bb3a8
Thu Jun 15 16:49:33 2017 +0200 ; vcfucsc ; https://github.com/lindenb/jvarkit/commit/66a41f0f51e057fa3b5c1281ff3d6539b12ae6f5
Thu Jun 15 15:30:26 2017 +0200 ; update vcfcalledwithanothermethod, vcfucsc ; https://github.com/lindenb/jvarkit/commit/0efbf47c1a7be8ee9b0a6e2e1dbfd82ae0f8508f
Sat Jun 3 23:36:42 2017 +0200 ; cleanup ; https://github.com/lindenb/jvarkit/commit/9574b7c9b25ef9d209f086f00e800481520cea67
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Fri Apr 7 16:35:31 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/54c5a476e62e021ad18e7fd0d84bf9e5396c8c96
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Mon Feb 3 18:12:01 2014 +0100 ; lundi. je rentre en velo ? il pleut... ; https://github.com/lindenb/jvarkit/commit/66c43aa46b61bbc7f037b1799be5871e82794ab2
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Sun Jan 5 16:10:56 2014 +0100 ; vcf set dict ; https://github.com/lindenb/jvarkit/commit/f023bc9b0685266627a260c67813e7b76d42bef1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfucsc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```
java -jar dist/vcfucsc.jar --table snp142 -e '${name}' input.vcf
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr3	124290753	.	G	C	579.77	.	UCSC_HG19_SNP142=rs145115089
chr3	124290943	.	A	G	491.77	.	UCSC_HG19_SNP142=rs7372055
chr3	124291069	.	G	A	266.77	.	UCSC_HG19_SNP142=rs7373767
chr3	124291171	.	C	CA	240.73	.	.
chr3	124291245	.	A	G	563.77	.	UCSC_HG19_SNP142=rs12695439
chr3	124291351	.	A	G	194.77	.	UCSC_HG19_SNP142=rs7613600
chr3	124291416	.	G	T	308.77	.	UCSC_HG19_SNP142=rs73189597
chr3	124291579	.	T	C	375.77	.	UCSC_HG19_SNP142=rs7649882
```


