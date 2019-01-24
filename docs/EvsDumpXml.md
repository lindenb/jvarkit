# EvsDumpXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Download data from EVS http://evs.gs.washington.edu/EVS as XML file.


## Usage

```
Usage: evsdumpxml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --out
       output filename. must end with.xml will be indexed with tribble
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -L
      limit to L records (for debugging)
      Default: -1
    -N
       download using a step of  'N' bases
      Default: 1000000
    -s
      sort data. (default if filename specified with -o)
      Default: false

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew evsdumpxml
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/evs2bed/EvsDumpXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/evs2bed/EvsDumpXml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **evsdumpxml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$  java -jar dist/evs2xml.jar -L 10000 -o test.xml
$ cat test.xml
``
```xml
<?xml version="1.0" encoding="UTF-8"?>
<evsData xmlns="http://webservice.evs.gs.washington.edu/">
  <snpList>
    <positionString>1:69428</positionString>
    <chrPosition>69428</chrPosition>
    <alleles>G/T</alleles>
    <uaAlleleCounts>G=313/T=6535</uaAlleleCounts>
    <aaAlleleCounts>G=14/T=3808</aaAlleleCounts>
    <totalAlleleCounts>G=327/T=10343</totalAlleleCounts>
    <uaMAF>4.5707</uaMAF>
    <aaMAF>0.3663</aaMAF>
    <totalMAF>3.0647</totalMAF>
    <avgSampleReadDepth>110</avgSampleReadDepth>
    <geneList>OR4F5</geneList>
    <snpFunction>
      <chromosome>1</chromosome>
      <position>69428</position>
      <conservationScore>1.0</conservationScore>
      <conservationScoreGERP>0.9</conservationScoreGERP>
      <snpFxnList>
        <mrnaAccession>NM_001005484.1</mrnaAccession>
        <fxnClassGVS>missense</fxnClassGVS>
        <hgvsProteinVar>p.(F113C)</hgvsProteinVar>
        <hgvsCdnaVar>c.338T&gt;G</hgvsCdnaVar>
        <codingDnaSize>918</codingDnaSize>
        <pphPrediction>probably-damaging:0.999</pphPrediction>
        <granthamScore>205</granthamScore>
      </snpFxnList>
      <refAllele>T</refAllele>
      <ancestralAllele>T</ancestralAllele>
      <firstRsId>140739101</firstRsId>
      <approxMapped2RsId>false</approxMapped2RsId>
      <filters>PASS</filters>
      <clinicalLink>unknown</clinicalLink>
    </snpFunction>
    <conservationScore>1.0</conservationScore>
    <conservationScoreGERP>0.9</conservationScoreGERP>
    <refAllele>T</refAllele>
    <altAlleles>G</altAlleles>
    <ancestralAllele>T</ancestralAllele>
    <chromosome>1</chromosome>
    <hasAtLeastOneAccession>true</hasAtLeastOneAccession>
    <rsIds>rs140739101</rsIds>
    <filters>PASS</filters>
    <clinicalLink>unknown</clinicalLink>
    <dbsnpVersion>dbSNP_134</dbsnpVersion>
    <uaGenotypeCounts>GG=92/GT=129/TT=3203</uaGenotypeCounts>
    <aaGenotypeCounts>GG=1/GT=12/TT=1898</aaGenotypeCounts>
    <totalGenotypeCounts>GG=93/GT=141/TT=5101</totalGenotypeCounts>
    <onExomeChip>false</onExomeChip>
    <gwasPubmedIds>unknown</gwasPubmedIds>
    <eaMutAge>-1.0</eaMutAge>
    <eaMutAgeSd>-1.0</eaMutAgeSd>
    <aaMutAge>-1.0</aaMutAge>
    <aaMutAgeSd>-1.0</aaMutAgeSd>
    <grcH38Position>1:69428</grcH38Position>
  </snpList>
  <snpList>
    <positionString>1:69476</positionString>
    <chrPosition>69476</chrPosition>
(...)
```
