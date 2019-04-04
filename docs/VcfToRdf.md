# VcfToRdf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert VCF to RDF (N3 notation)


## Usage

```
Usage: vcf2rdf [options] Files
  Options:
    -a, --alleles
      print ALT alleles
      Default: false
    -f, --filters
      print FILTERs
      Default: false
    -g, --genotypes
      print Genotypes informations
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -vep, --vep
      print VEP informations
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * rdf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcf2rdf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2rdf/VcfToRdf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2rdf/VcfToRdf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2rdf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$  java -jar dist/vcf2rdf.jar < in.vcf | xmllint --format -
```


```xml
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:vcf="http://github.com/lindenb/jvarkit/" xmlns:xsd="http://www.w3.org/2001/XMLSchema#">
  <vcf:Chromosome rdf:about="urn:chromosome/1">
    <dc:title>1</dc:title>
    <vcf:length rdf:datatype="xsd:int">249250621</vcf:length>
    <vcf:index rdf:datatype="xsd:int">0</vcf:index>
  </vcf:Chromosome>
  <vcf:Chromosome rdf:about="urn:chromosome/2">
    <dc:title>2</dc:title>
    <vcf:length rdf:datatype="xsd:int">243199373</vcf:length>
(...)
  <vcf:Chromosome rdf:about="urn:chromosome/GL000192.1">
    <dc:title>GL000192.1</dc:title>
    <vcf:length rdf:datatype="xsd:int">547496</vcf:length>
    <vcf:index rdf:datatype="xsd:int">83</vcf:index>
  </vcf:Chromosome>
  <vcf:Filter rdf:about="urn:filter/FILTER">
    <dc:title>FILTER</dc:title>
    <dc:description/>
  </vcf:Filter>
  <vcf:Sample rdf:about="urn:sample/V2528">
    <dc:title>V2528</dc:title>
  </vcf:Sample>
  <vcf:Variant rdf:about="urn:variant/1">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467571</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467571</vcf:end>
    <vcf:ref>C</vcf:ref>
    <vcf:alt>C</vcf:alt>
    <vcf:alt>T</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">1311.77</vcf:qual>
    <vcf:BaseQRankSum>-2.441</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>365.0758</vcf:HaplotypeScore>
    <vcf:QD>5.33</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>49.60</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>88.037</vcf:FS>
    <vcf:MQRankSum>-10.623</vcf:MQRankSum>
    <vcf:ReadPosRankSum>-5.556</vcf:ReadPosRankSum>
    <vcf:Dels>0.02</vcf:Dels>
    <vcf:DP>246</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/2">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/1"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>C</vcf:allele>
    <vcf:allele>T</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">229</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>1340,0,5602</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/3">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467576</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467576</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">624.77</vcf:qual>
    <vcf:BaseQRankSum>-4.704</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>415.3829</vcf:HaplotypeScore>
    <vcf:QD>2.58</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>48.90</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>26.182</vcf:FS>
    <vcf:MQRankSum>4.246</vcf:MQRankSum>
    <vcf:ReadPosRankSum>0.781</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>242</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/4">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/3"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">230</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>653,0,6159</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/5">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467578</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467581</vcf:end>
    <vcf:ref>TTTC</vcf:ref>
    <vcf:alt>TTTC</vcf:alt>
    <vcf:alt>T</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">2.14748360973E9</vcf:qual>
    <vcf:FS>7.426</vcf:FS>
    <vcf:AC>1</vcf:AC>
    <vcf:BaseQRankSum>-3.400</vcf:BaseQRankSum>
    <vcf:MQRankSum>-7.221</vcf:MQRankSum>
    <vcf:ReadPosRankSum>-3.139</vcf:ReadPosRankSum>
    <vcf:DP>243</vcf:DP>
    <vcf:QD>29.69</vcf:QD>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>49.83</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/6">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/5"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>TTTC</vcf:allele>
    <vcf:allele>T</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">243</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
   <h:pre><![CDATA[ <vcf:pl>5687,0,7242</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/7">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467588</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467588</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">610.77</vcf:qual>
    <vcf:BaseQRankSum>5.163</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>106.1491</vcf:HaplotypeScore>
    <vcf:QD>2.51</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>54.46</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>18.558</vcf:FS>
    <vcf:MQRankSum>3.886</vcf:MQRankSum>
    <vcf:ReadPosRankSum>0.762</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>243</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/8">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/7"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">231</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>639,0,6257</vcf:pl>
  </vcf:Genotype>
  <vcf:Variant rdf:about="urn:variant/9">
    <vcf:chromosome rdf:resource="urn:chromosome/8"/>
    <vcf:start rdf:datatype="xsd:int">10467589</vcf:start>
    <vcf:end rdf:datatype="xsd:int">10467589</vcf:end>
    <vcf:ref>T</vcf:ref>
    <vcf:alt>T</vcf:alt>
    <vcf:alt>C</vcf:alt>
    <vcf:qual rdf:datatype="xsd:double">3705.77</vcf:qual>
    <vcf:BaseQRankSum>6.528</vcf:BaseQRankSum>
    <vcf:HaplotypeScore>79.8219</vcf:HaplotypeScore>
    <vcf:QD>15.19</vcf:QD>
    <vcf:MLEAC>1</vcf:MLEAC>
    <vcf:MQ>55.35</vcf:MQ>
    <vcf:MLEAF>0.500</vcf:MLEAF>
    <vcf:AC>1</vcf:AC>
    <vcf:FS>19.873</vcf:FS>
    <vcf:MQRankSum>2.222</vcf:MQRankSum>
    <vcf:ReadPosRankSum>4.388</vcf:ReadPosRankSum>
    <vcf:Dels>0.00</vcf:Dels>
    <vcf:DP>244</vcf:DP>
    <vcf:AF>0.500</vcf:AF>
    <vcf:MQ0>0</vcf:MQ0>
    <vcf:AN>2</vcf:AN>
  </vcf:Variant>
  <vcf:Genotype rdf:about="urn:genotype/10">
    <vcf:sample rdf:resource="urn:sample/V2528"/>
    <vcf:variant rdf:resource="urn:variant/9"/>
    <rdf:type rdf:resource="urn:genotype/het"/>
    <vcf:allele>T</vcf:allele>
    <vcf:allele>C</vcf:allele>
    <vcf:dp rdf:datatype="xsd:int">232</vcf:dp>
    <vcf:gq rdf:datatype="xsd:int">99</vcf:gq>
    <vcf:pl>3734,0,3449</vcf:pl>
  </vcf:Genotype>
</rdf:RDF>
```


### Example

load with jena TDB


```

(cd /home/lindenb/src/jvarkit-git/ && make vcf2rdf)
export JENAROOT=/home/lindenb/package/apache-jena-2.11.0
 
rm -rf TMPTDB
for F in *.vcf.gz
do
gunzip -c ${F} |\
java -jar ${HOME}/src/jvarkit-git/dist/vcf2rdf.jar > tmp.rdf
${JENAROOT}/bin/tdbloader2 --loc=TMPTDB tmp.rdf
du -hs TMPTDB
done
rm tmp.rdf
${JENAROOT}/bin/tdbdump --loc=TMPTDB | head -n 100

loaded 78,258 variants / 20 Mbytes in Jena/TDB 11,193,250 triples / 892Mbytes 

```





