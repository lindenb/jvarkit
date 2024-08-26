# DictToXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert a SAM dictionary from vcf,sam,bam,dict, etc.. to XML.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar dict2xml  [options] Files

Usage: dict2xml [options] Files
  Options:
    --duplicate
      keep duplicates (default behavior is to keep one dictionary
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-errors
      ignore errors, skip files that don't have a dictionary
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * dict
 * xml
 * sam
 * bam
 * vcf



## Creation Date

20240824

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2xml/DictToXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2xml/DictToXml.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **dict2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to XML
Then we can use XSLT to generate code...

## Example

```
$ java -jar dist/jvarkit.jar dict2xml "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" | xmllint --format - | head
<?xml version="1.0" encoding="UTF-8"?>
<dictionaries>
  <dictionary md5="c4e11bf85a6ab2f944340f409c751f32" length="3137454505" count="86" source="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
    <sequence name="2" length="243199373" index="1" offset="249250621"/>
    <sequence name="3" length="198022430" index="2" offset="492449994"/>
    <sequence name="4" length="191154276" index="3" offset="690472424"/>
    <sequence name="5" length="180915260" index="4" offset="881626700"/>
    <sequence name="6" length="171115067" index="5" offset="1062541960"/>
    <sequence name="7" length="159138663" index="6" offset="1233657027"/>
(...)
```

```
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2xml | xmllint --format - 
<?xml version="1.0" encoding="UTF-8"?>
<dictionaries>
  <dictionary md5="4677ece43eea2b029d0d33fe130ea6c7" length="3137454505" count="86" source="src/test/resources/roxan
.hs37d5.csq.vcf.gz" build="GRCh37">
    <sequence name="chr1" length="249250621" index="0" offset="0"/>
(...)
    <sequence name="hs37d5" length="35477943" index="85" offset="3101976562"/>
  </dictionary>
  <dictionary md5="bd7e0928fc3c810e48fafc53a4222ed5" length="18490" count="11" source="src/test/resources/S4.vcf.gz"
>
    <sequence name="RF01" length="3302" index="0" offset="0"/>
(...)
    <sequence name="RF10" length="751" index="9" offset="17073"/>
    <sequence name="RF11" length="666" index="10" offset="17824"/>
  </dictionary>
  <dictionary md5="9a5c58c2c91e731135b27ed14974523a" length="3101976562" count="85" source="src/test/resources/gnoma
d.exomes.r2.0.1.sites.vcf.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
(...)
    <sequence name="NC_007605" length="171823" index="84" offset="3101804739"/>
  </dictionary>
  <dictionary md5="635de5cb51973d45844fa713ac0b7719" length="85" count="2" source="src/test/resources/toy.vcf.gz">
    <sequence name="ref" length="45" index="0" offset="0"/>
    <sequence name="ref2" length="40" index="1" offset="45"/>
  </dictionary>
  <dictionary md5="f8d942cb3fc6ebef618a0b0ba3f4ef99" length="3095677412" count="24" source="src/test/resources/gnoma
d_v2_sv.sites.vcf.gz" build="GRCh37">
    <sequence name="1" length="249250621" index="0" offset="0"/>
    <sequence name="2" length="243199373" index="1" offset="249250621"/>
    <sequence name="3" length="198022430" index="2" offset="492449994"/>
    <sequence name="4" length="191154276" index="3" offset="690472424"/>
(...)
    <sequence name="Y" length="59373566" index="23" offset="3036303846"/>
  </dictionary>
</dictionaries>

```


