# VCFBigWig

Annotate a VCF with values from a bigwig file


## Usage

```
Usage: vcfbigwig [options] Files
  Options:
    -a, --aggregate
      How to aggregate overlapping values: 'avg' average; 'median': median, 
      'first': use first, 'all' : print all the data
      Default: avg
      Possible Values: [avg, median, first, all]
  * -B, --bigwig
      Path to the bigwig file
    -C, --contained
      Specifies wig values must be contained by region. if false: return any 
      intersecting region values
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --onNotFound
      [20170707] Contig converter. I will do my best to convert the contig 
      names (e.g 'chr1' -> '1'): But what should I do when comparing two 
      dictionaries with different notations
      Default: SKIP
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tag, -tag
      Name of the INFO tag. default: name of the bigwig
    --version
      print version and exit

```


## Keywords

 * vcf
 * wig
 * wiggle
 * bigwig


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
$ make vcfbigwig
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWig.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbigwig/VCFBigWig.java)

Git History for this file:
```
Fri Sep 8 11:29:13 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/0c5035b59124ca18aba0405d0c616b565a32d10e
Wed Sep 6 14:49:24 2017 +0200 ; fixing typos, starting to generate VariantContextWriterFactory for spring xml ; https://github.com/lindenb/jvarkit/commit/cf023e059af85f6c266c56a8f7db6ff78e4a5134
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Thu Jul 27 16:58:18 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/a8aaf2d7df89f44442b36ee1120ee4dd5c1e36e6
Fri Jul 7 18:36:14 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c5dc2be25578f7cbc60c0f5425bacf4450893c92
Mon May 29 12:33:45 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/870be8e90d7e98d947f73e67ef9965f12f351846
Sun May 21 20:02:10 2017 +0200 ; instanceMain -> instanceMainWithExit ; https://github.com/lindenb/jvarkit/commit/4fa41d198fe7e063c92bdedc333cbcdd2b8240aa
Fri May 19 17:10:13 2017 +0200 ; cont doc ; https://github.com/lindenb/jvarkit/commit/d2aea1eaa554d0498b197fb8fac01893b10ceb83
Fri Apr 21 18:16:07 2017 +0200 ; scan sv ; https://github.com/lindenb/jvarkit/commit/49b99018811ea6a624e3df556627ebdbf3f16eab
Fri Mar 17 17:06:05 2017 +0100 ; added SO terms, added vep to bioalcidae ; https://github.com/lindenb/jvarkit/commit/96d5c5dcc556f399bb8cf34bbc4d6f31fbebc8c5
Wed Mar 15 17:49:49 2017 +0100 ; fix window, enhance vcfbigwig ; https://github.com/lindenb/jvarkit/commit/af350cdad1b64edfcd2b02f1830b814c9115e31a
Mon Jan 18 16:58:08 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/83f80fdbe8d6be71539cfdbf60d61ce7ead9c0fd
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Thu Apr 30 16:54:43 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/eb6a52869eb3b9b8bf048a079e1d7a96bccab2bf
Tue Apr 28 15:22:51 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/10924ba7a0d78e208109157b92432d6171247c21
Fri Apr 24 17:05:31 2015 +0200 ; moved VCF-BigWig to a standard command-line #tweet ; https://github.com/lindenb/jvarkit/commit/3a0c4ccb05e7492382e00328ac60951f215d9400
Tue Feb 24 16:43:03 2015 +0100 ; vcfin : code rewrittern. picky with ALT alleles. #tweet ; https://github.com/lindenb/jvarkit/commit/65ef7741539e89c7a1a1f9cca28c13d531902c96
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Mon Jan 13 12:47:15 2014 +0100 ; vcf bigwig ; https://github.com/lindenb/jvarkit/commit/ffcdd7bf241c12030cc6b67fa912429a0c9b5e31
Mon Jan 13 11:19:02 2014 +0100 ; Making contained a option ; https://github.com/lindenb/jvarkit/commit/b27f41e37d4a0cb5274081ec9d850ee2ac5bf014
Fri Oct 11 15:39:02 2013 +0200 ; picard v.100: deletion of VcfIterator :-( ; https://github.com/lindenb/jvarkit/commit/e88fab449b04aed40c2ff7f9d0cf8c8b6ab14a31
Fri Sep 6 15:11:11 2013 +0200 ; moved code for latest version of picard (1.97). Using VCFIterator instead of ASciiLineReader ; https://github.com/lindenb/jvarkit/commit/810877c10406a017fd5a31dacff7e8401089d429
Fri Jul 26 10:38:30 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/1d4be124f6de986aa54102970d1ea1a32d8f8329
Wed Jul 24 14:44:13 2013 +0200 ; starting my code from http://plindenbaum.blogspot.com/2011/01/my-tool-to-annotate-vcf-files.html ; https://github.com/lindenb/jvarkit/commit/cbd77f6fc09edd992990112cbfc959b3b09574d5
Sun Jul 21 14:17:59 2013 +0200 ; vcf trios, added git HASH in METAINF/Manifest ; https://github.com/lindenb/jvarkit/commit/1854d3695563b91471861164f5e8903042493470
Wed Jul 17 08:39:12 2013 +0200 ; tmp ; https://github.com/lindenb/jvarkit/commit/8dc0f2623a4e6e5a5c849ca8083efe689819fea9
Tue Jul 16 13:13:34 2013 +0200 ; moving bigwig to picard ; https://github.com/lindenb/jvarkit/commit/6b32bcf0385fa6d1125b97a7722cf99c82f7ead4
Mon May 6 18:56:46 2013 +0200 ; moving to git ; https://github.com/lindenb/jvarkit/commit/55158d13f0950f16c4a3cc3edb92a87905346ee1
```

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbigwig** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
 java -jar dist/vcfbigwig.jar \
 	-T GERP \
 	-B gerp.bw input.vcf.gz 
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```



