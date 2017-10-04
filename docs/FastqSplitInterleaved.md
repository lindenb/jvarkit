# FastqSplitInterleaved

Split interleaved Fastq files.


## Usage

```
Usage: fastqsplitinterleaved [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --version
      print version and exit
    -a
      (fastq1 file or '-' for stdout). Ignore 1st read if omitted. Optional.
    -b
      (fastq2 file or '-' for stdout). Ignore 2nd read if omitted. Optional.

```


## Keywords

 * fastq


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
$ make fastqsplitinterleaved
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplitInterleaved.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqSplitInterleaved.java)


<details>
<summary>Git History</summary>

```
Fri May 12 19:41:30 2017 +0200 ; fix make, empty doc ; https://github.com/lindenb/jvarkit/commit/52fcf6d46a779fd7153ebc032fae643d2e266e7e
Fri May 12 18:07:46 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ca96bce803826964a65de33455e5231ffa6ea9bd
Tue Apr 18 13:24:50 2017 +0200 ; cont-cleanup ; https://github.com/lindenb/jvarkit/commit/a86c8971fe5ebb3f8de175c75e78f2d0e5325cfd
Fri Jun 17 13:56:39 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/865252a44fc018f46b4280788cec65a1383dcc18
Mon Sep 1 15:47:33 2014 +0200 ; fix error in shuffle. pubmed filter js ; https://github.com/lindenb/jvarkit/commit/1c690742a61ba809e342eacd5d2a214134bfab72
Mon Sep 1 12:18:15 2014 +0200 ; split interleaved fastq ; https://github.com/lindenb/jvarkit/commit/5be82f90e563427ea157e3ecdbba1922a0dc37a3
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqsplitinterleaved** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ curl -sk "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
java -jar dist/fastqsplitinterleaved.jar  -a -  | grep "^@"

@H06HDADXX130110:2:2116:3345:91806/1
@H06HDADXX130110:1:2103:11970:57672/1
@H06JUADXX130110:1:1108:6424:55322/1

$ curl -sk "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
java -jar dist/fastqsplitinterleaved.jar  -b -  | grep "^@"
@H06HDADXX130110:2:2116:3345:91806/2
@H06HDADXX130110:1:2103:11970:57672/2
@H06JUADXX130110:1:1108:6424:55322/2

$ curl -sk "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
java -jar dist/fastqsplitinterleaved.jar  -a x1.fq.gz -b x2.fq.gz


