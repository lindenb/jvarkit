# Biostar59647

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

SAM/BAM to XML


## Usage

```
Usage: biostar59647 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
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

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar59647
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar59647Test.java)


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

