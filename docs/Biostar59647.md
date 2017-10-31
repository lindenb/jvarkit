# Biostar59647

SAM/BAM to XML


## Usage

```
Usage: biostar59647 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/59647](https://www.biostars.org/p/59647)


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) and avoid OpenJdk, use the java from Oracle. Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make biostar59647
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647.java)


<details>
<summary>Git History</summary>

```
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Fri Apr 14 15:27:32 2017 +0200 ; annotation proc ; https://github.com/lindenb/jvarkit/commit/72b9383a8472e5a91120bab84d15b8acad4db8d4
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Sat Dec 14 00:39:57 2013 +0100 ; updated Biostar59647 ; https://github.com/lindenb/jvarkit/commit/35f92334d89c1e970b1e9c0f9f076ade021d5492
Tue Nov 12 13:18:49 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/7fa73b9b3d540c6fe444517ac10bf4ee945a03e1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar59647** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ java -jar dist/biostar59647.jar -r samtools-0.1.18/examples/toy.fa  samtools-0.1.18/examples/toy.bam |\
xmllint --format - 
```
```xml
<?xml version="1.0" encoding="UTF-8"?>
<sam ref="/home/lindenb/samtools-0.1.18/examples/toy.bam" bam="/home/lindenb/samtools-0.1.18/examples/toy.fa">
  <!---r /home/lindenb/samtools-0.1.18/examples/toy.fa /home/lindenb/samtools-0.1.18/examples/toy.bam-->
  <read>
    <name>r001</name>
    <sequence>TTAGATAAAGAGGATACTG</sequence>
    <flags READ_PAIRED="true" READ_MAPPED_IN_PROPER_PAIR="true" READ_UNMAPPED="false" MATE_UNMAPPED="false" READ
_REVERSE_STRAND="false" MATE_REVERSE_STRAND="true" FIRST_IN_PAIR="false" SECOND_IN_PAIR="true" NOT_PRIMARY_ALIGN
MENT="false" READ_FAILS_VENDOR_QUALITY_CHECK="false" READ_IS_DUPLICATE="false" SUPPLEMENTARY_ALIGNMENT="false">1
63</flags>
    <qual>30</qual>
    <chrom index="0">ref</chrom>
    <pos>7</pos>
    <cigar>8M4I4M1D3M</cigar>
    <mate-chrom index="0">ref</mate-chrom>
    <mate-pos>37</mate-pos>
    <align>
      <M read-index="1" read-base="T" ref-index="7" ref-base="T"/>
      <M read-index="2" read-base="T" ref-index="8" ref-base="T"/>
      <M read-index="3" read-base="A" ref-index="9" ref-base="A"/>
      <M read-index="4" read-base="G" ref-index="10" ref-base="G"/>
      <M read-index="5" read-base="A" ref-index="11" ref-base="A"/>
      <M read-index="6" read-base="T" ref-index="12" ref-base="T"/>
      <M read-index="7" read-base="A" ref-index="13" ref-base="A"/>
      <M read-index="8" read-base="A" ref-index="14" ref-base="A"/>
      <I read-index="9" read-base="A"/>
      <I read-index="10" read-base="G"/>
      <I read-index="11" read-base="A"/>
      <I read-index="12" read-base="G"/>
      <M read-index="13" read-base="G" ref-index="15" ref-base="G"/>
      <M read-index="14" read-base="A" ref-index="16" ref-base="A"/>
      <M read-index="15" read-base="T" ref-index="17" ref-base="T"/>
      <M read-index="16" read-base="A" ref-index="18" ref-base="A"/>
      <D ref-index="19" ref-base="G"/>
      <M read-index="17" read-base="C" ref-index="20" ref-base="C"/>
      <M read-index="18" read-base="T" ref-index="21" ref-base="T"/>
      <M read-index="19" read-base="G" ref-index="22" ref-base="G"/>
    </align>
  </read>
  <read>
    <name>r002</name>
    <sequence>AAAAGATAAGGGATAAA</sequence>
    <flags READ_PAIRED="false" READ_MAPPED_IN_PROPER_PAIR="false" READ_UNMAPPED="false" MATE_UNMAPPED="false" RE
AD_REVERSE_STRAND="false" MATE_REVERSE_STRAND="false" FIRST_IN_PAIR="false" SECOND_IN_PAIR="false" NOT_PRIMARY_A
LIGNMENT="false" READ_FAILS_VENDOR_QUALITY_CHECK="false" READ_IS_DUPLICATE="false" SUPPLEMENTARY_ALIGNMENT="fals
e">0</flags>
    <qual>30</qual>
    <chrom index="0">ref</chrom>
    <pos>9</pos>
    <cigar>1S2I6M1P1I1P1I4M2I</cigar>
    <align>
      <S read-index="1" read-base="A"/>
      <I read-index="2" read-base="A"/>
      <I read-index="3" read-base="A"/>
      <M read-index="4" read-base="A" ref-index="9" ref-base="A"/>
      <M read-index="5" read-base="G" ref-index="10" ref-base="G"/>
      <M read-index="6" read-base="A" ref-index="11" ref-base="A"/>
      <M read-index="7" read-base="T" ref-index="12" ref-base="T"/>
      <M read-index="8" read-base="A" ref-index="13" ref-base="A"/>
      <M read-index="9" read-base="A" ref-index="14" ref-base="A"/>
      <I read-index="10" read-base="G"/>
      <I read-index="11" read-base="G"/>
      <M read-index="12" read-base="G" ref-index="15" ref-base="G"/>
(...)
```

## Cited in

* cited in http://biorxiv.org/content/early/2014/01/21/001834 "Illumina TruSeq synthetic long-reads empower de novo assembly and resolve complex, highly repetitive transposable elements"


