# GroupByGene

Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff


## Usage

```
Usage: groupbygene [options] Files
  Options:
    --annIntergenic
      [20170726] Accept snpEff 'ANN' intergenic regions.
      Default: false
    --filtered
      ignore FILTERED variants
      Default: false
    --fisher
      [20170726] Print fisher for case/control (experimental, need to work on 
      this) 
      Default: false
    --gtFiltered
      [20170725] ignore FILTERED genotypes
      Default: false
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
    -p, --ped, --pedigree
      [20170725] A pedigree is a text file delimited with tabs. No header. 
      Columns are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother 
      Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 
      unaffected, 1 affected,-9 unknown
    --slidingWindowShift
      [20170726] if greater than 0, add shift the sliding window by this 
      distance 
      Default: 0
    --slidingWindowSize
      [20170726] if greater than 0, add a sliding window of this size as the 
      name of a virtual region
      Default: 0
    -T, --tag
      search Tag in INFO field containing the name of the genes (optional)
      Default: []
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --typeRegexExclude
      [20170726] ignore prediction type matching this regex. e.g: 
      '(ann_gene_id|ann_feature_id)' 
    --version
      print version and exit

```


## Keywords

 * vcf
 * gene


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
$ make groupbygene
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java)


<details>
<summary>Git History</summary>

```
Mon Aug 7 09:53:19 2017 +0200 ; fixed unicode problems after https://github.com/lindenb/jvarkit/issues/82 ; https://github.com/lindenb/jvarkit/commit/68254c69b027a9ce81d8b211447f1c0bf02dc626
Wed Jul 26 18:09:38 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/576fdd17812f9a47491945cb8bb74990ffb084c9
Tue Jul 25 11:32:47 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7f9befb99c8771e45e437e7d7f178195d15fdbaa
Wed Jun 28 17:33:30 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/3c252f19e5cad0ec87d250a5b9884b6f2d6fe856
Mon Jun 26 17:29:03 2017 +0200 ; burden ; https://github.com/lindenb/jvarkit/commit/a3b7abf21d07f0366e81816ebbb2cce26b2341e7
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Wed Apr 5 13:49:50 2017 +0200 ; cont, fix bug in findallcovatpos ; https://github.com/lindenb/jvarkit/commit/7db18c7fe90fd5bf64d3ff3a4505607a1974ce6b
Tue Apr 4 17:09:36 2017 +0200 ; vcfgnomad ; https://github.com/lindenb/jvarkit/commit/eac33a01731eaffbdc401ec5fd917fe345b4a181
Wed Feb 22 19:07:03 2017 +0100 ; refactor prediction parsers ; https://github.com/lindenb/jvarkit/commit/dc7f7797c60d63cd09d3b7712fb81033cd7022cb
Wed Feb 3 14:41:34 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/dea087858303eb791d6e68178742f1fbae2092f0
Fri Jan 29 10:46:18 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/f94556b6ac4d0f0828c83b5b6bcb9cdab746979e
Fri Jan 22 23:53:30 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/88b87891d1f2336c4a76c34ae34728203aa8ac0f
Fri Jan 22 23:49:23 2016 +0100 ; vcfiterator is now an interface ; https://github.com/lindenb/jvarkit/commit/9f9b9314c4b31b21044c5911a7e79e1b3fb0af7a
Mon Nov 30 16:53:51 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/89f3cbe043ac8c52735feec5b45e43cf873b7179
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Thu Feb 26 17:31:31 2015 +0100 ; automatic generation of tool.xml for #usegalaxy using xslt. #tweet ; https://github.com/lindenb/jvarkit/commit/d3074f97acbc2ccf095efcc587c6faf7cffa7d03
Mon Jan 26 13:02:23 2015 +0100 ; (tweet) updated groupbygene with latest version of VEP ; https://github.com/lindenb/jvarkit/commit/6dd9f22bf5a700dfaaa044b9ff3e86590583beb5
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Fri Dec 13 17:19:27 2013 +0100 ; fix major bug in prediction parsers ; https://github.com/lindenb/jvarkit/commit/7559f47ef5f1ecee018ea7eb0968b6bdede93283
Tue Dec 10 12:30:02 2013 +0100 ; cont. ; https://github.com/lindenb/jvarkit/commit/88a74712aabfc418d7ed390fff6e965e7491b210
Mon Dec 9 14:47:17 2013 +0100 ; group by gene ; https://github.com/lindenb/jvarkit/commit/5f32d705f9e553cde494610def4e6c7f98dc5412
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **groupbygene** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Group VCF data by gene/transcript. By default it tries to use data from VEP and SnpEff

## Example

### Delimited output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar |\
head | column  -t

#chrom  min.POS    max.POS    gene.name  gene.type         samples.affected  count.variations  M10475  M10478  M10500  M128215
chr10   52004315   52004315   ASAH2      snpeff-gene-name  2                 1                 0       0       1       1
chr10   52004315   52004315   ASAH2      vep-gene-name     2                 1                 0       0       1       1
chr10   52497529   52497529   ASAH2B     snpeff-gene-name  2                 1                 0       1       1       0
chr10   52497529   52497529   ASAH2B     vep-gene-name     2                 1                 0       1       1       0
chr10   48003992   48003992   ASAH2C     snpeff-gene-name  3                 1                 1       1       1       0
chr10   48003992   48003992   ASAH2C     vep-gene-name     3                 1                 1       1       1       0
chr10   126678092  126678092  CTBP2      snpeff-gene-name  1                 1                 0       0       0       1
chr10   126678092  126678092  CTBP2      vep-gene-name     1                 1                 0       0       0       1
chr10   135336656  135369532  CYP2E1     snpeff-gene-name  3                 2                 0       2       1       1
```

```

## History

* 201707: added pedigree, removed XML output



