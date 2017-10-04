# GcPercentAndDepth

Extracts GC% and depth for multiple bam using a sliding window


## Usage

```
Usage: gcpercentanddepth [options] Files
  Options:
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
    --version
      print version and exit
    -B
       (file.bed) (optional). If not defined: use whole genome. Warning memory 
      consumming: must alloc sizeof(int)*win.size()*num(samples).
    -N
       (file) . chrom.name.helper .
    -R
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -m
       min depth
      Default: 0
    -n
       skip window if Reference contains one 'N'.
      Default: false
    -o
      Output file. Optional . Default: stdout
    -s
       (window shift)
      Default: 50
    -w
       (window size)
      Default: 100
    -x
       don't print genomic index.
      Default: false

```


## Keywords

 * gc%
 * depth
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
$ make gcpercentanddepth
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/GcPercentAndDepth.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/GcPercentAndDepth.java)


<details>
<summary>Git History</summary>

```
Fri Jun 2 16:31:30 2017 +0200 ; circos / lumpy ; https://github.com/lindenb/jvarkit/commit/7bddffca3899196e568fb5e1a479300c0038f74f
Mon May 29 12:33:45 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/870be8e90d7e98d947f73e67ef9965f12f351846
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Mon May 15 10:41:51 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c13a658b2ed3bc5dd6ade57190e1dab05bf70612
Thu May 4 13:06:07 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/b2f8f945cb8838c0289a7d850ce24603417eccde
Wed Apr 22 12:21:29 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/dc65752a1d0c364957940847f8901d32106f21c7
Mon Apr 13 10:32:12 2015 +0200 ; added option -n for https://github.com/lindenb/jvarkit/issues/26#issuecomment-92193999 ; https://github.com/lindenb/jvarkit/commit/19072426981ab0f13e61755012a4a35da591cc95
Tue Feb 24 16:43:03 2015 +0100 ; vcfin : code rewrittern. picky with ALT alleles. #tweet ; https://github.com/lindenb/jvarkit/commit/65ef7741539e89c7a1a1f9cca28c13d531902c96
Thu Nov 27 13:11:06 2014 +0100 ; bam compare coverage ; https://github.com/lindenb/jvarkit/commit/0be60cca2b40fa2bb2713e759271573936911aba
Thu Jun 5 11:06:35 2014 +0200 ; messages ; https://github.com/lindenb/jvarkit/commit/a5d02dbfcb811e9ed577a8850033c86ef5a91332
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Apr 9 12:50:50 2014 +0200 ; bwa mem convert cigar to N ; https://github.com/lindenb/jvarkit/commit/dac74ffe48ee2fbf391a9e8e0056e9b8ee18babc
Wed Apr 2 17:50:50 2014 +0200 ; cont rnaseq ; https://github.com/lindenb/jvarkit/commit/7b3f7e13a112b09018284931678ac78dd32cefcc
Wed Mar 5 15:25:14 2014 +0100 ; continue everything ; https://github.com/lindenb/jvarkit/commit/15fc86fe1703398069d0bba3300be31c396dc818
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gcpercentanddepth** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$ java -jar dist/gcanddepth.jar -R ref.fasta -b capture.bed 1.bam 2.bam ... > result.tsv
```

