# GroupByGene


## Usage

```
Usage: groupbygene [options] Files
  Options:
    ---maxRecordsInRam
      Max records in RAM
      Default: 50000
    --filtered
      ignore FILTERED variants
      Default: false
    -h, --help
      print help and exit
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tag
      add Tag in INFO field containing the name of the genes.
      Default: []
    --tmpDir
      Temporary directory
      Default: /tmp
    --version
      print version and exit
    -X, --xml
      XML output
      Default: false

```


## Description

Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff


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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java
)
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

### Delimited output

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

## XML output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar -X |\
xmllint --> --format -<!-- 
```


```
<!-- <?xml version="1.0" encoding="UTF-8"?>
<genes>
  <samples count="4">
    <sample>M10475</sample>
    <sample>M10478</sample>
    <sample>M10500</sample>
    <sample>M128215</sample>
  </samples>
  <gene name="ASAH2" type="snpeff-gene-name" chrom="chr10" min.POS="52004315" max.POS="52004315" affected="2" variations="1">
    <sample name="M10500" count="1">
      <genotype pos="52004315" ref="T" A1="C" A2="C"/>
    </sample>
    <sample name="M128215" count="1">
      <genotype pos="52004315" ref="T" A1="C" A2="C"/>
    </sample>
  </gene>
  <gene name="ASAH2" type="vep-gene-name" chrom="chr10" min.POS="52004315" max.POS="52004315" affected="2" variations="1">
    <sample name="M10500" count="1">
(...)
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000572003" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000572887" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000573843" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000573922" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000574309" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
</genes>
```


