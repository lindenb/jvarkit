# VcfSamplesPRS

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

another program for @AntoineRimbert


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfsamplesprs  [options] Files

Usage: vcfsamplesprs [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * --score, --scores, -S
      tabix indexed scores. 
      CHROM(tab)POS(tab)ID(tab)REF(tab)ALT(tab)EFFECT_ALL(tab)EFFECT 
    --version
      print version and exit

```


## Keywords

 * vcf
 * indel



## Creation Date

20220915

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/prs/VcfSamplesPRS.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/prs/VcfSamplesPRS.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsamplesprs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# buid tabix database

```
$ head AR_test.txt | column -t
#CHROM  POS        ID          REF  ALT  EFFECT_ALL  EFFECT  
chr1    5550460   rs249409   G    A    G           0.052
chr1    10918306  rs69301    G    T    T           0.15
chr2    2126390   rs13617   G    A    A           0.1
chr2    4407256   rs476   G    T    G           0.071
chr6    1605860  rs1564348   T    C    T           0.014
chr6    2093141   rs18   G    A    G           0.057

LC_ALL=C sort -T . -t $'\t' -k1,1 -k2,2n  AR_test.txt | bgzip > AR_test.txt.gz
tabix -s1 -b 2 -e 2  -c '#' 20220915_AR_test.txt.gz
```

and then

```
bcftools view in.bcf | java -jar dist/vcfsamplesprs.jar -S AR_test.txt.gz
```


