# SamFixCigar

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fix Cigar String in SAM replacing 'M' by 'X' or '='


## Usage

```
Usage: samfixcigar [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * cigar



## See also in Biostars

 * [https://www.biostars.org/p/312430](https://www.biostars.org/p/312430)
 * [https://www.biostars.org/p/340479](https://www.biostars.org/p/340479)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samfixcigar
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samfixcigar/SamFixCigar.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samfixcigar/SamFixCigar.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samfixcigar/SamFixCigarTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samfixcigar/SamFixCigarTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samfixcigar** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


the input file


```
$ cat toy.sam

@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r002    0       ref     9       30      1S2I6M1P1I1P1I4M2I      *       0       0       AAAAGATAAGGGATAAA       *
r003    0       ref     9       30      5H6M    *       0       0       AGCTAA  *
r004    0       ref     16      30      6M14N1I5M       *       0       0       ATAGCTCTCAGC    *
r003    16      ref     29      30      6H5M    *       0       0       TAGGC   *
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *
x1      0       ref2    1       30      20M     *       0       0       aggttttataaaacaaataa    ????????????????????
x2      0       ref2    2       30      21M     *       0       0       ggttttataaaacaaataatt   ?????????????????????
x3      0       ref2    6       30      9M4I13M *       0       0       ttataaaacAAATaattaagtctaca      ??????????????????????????
x4      0       ref2    10      30      25M     *       0       0       CaaaTaattaagtctacagagcaac       ?????????????????????????
x5      0       ref2    12      30      24M     *       0       0       aaTaattaagtctacagagcaact        ????????????????????????
x6      0       ref2    14      30      23M     *       0       0       Taattaagtctacagagcaacta ???????????????????????

```


processing with samfixcigar


```
$ java -jar dist/samfixcigar.jar \
     -r samtools-0.1.19/examples/toy.fa \
     samtools-0.1.19/examples/toy.sam
@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
r001    163     ref     7       30      8=4I4=1D3=      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r002    0       ref     9       30      1S2I6=1P1I1P1I1X1=2X2I  *       0       0       AAAAGATAAGGGATAAA       *
r003    0       ref     9       30      2=1X3=  *       0       0       AGCTAA  *
r004    0       ref     16      30      6=14N1I5=       *       0       0       ATAGCTCTCAGC    *
r003    16      ref     29      30      5=      *       0       0       TAGGC   *
r001    83      ref     37      30      9=      =       7       -39     CAGCGCCAT       *
x1      0       ref2    1       30      16=1X3= *       0       0       AGGTTTTATAAAACAAATAA    ????????????????????
x2      0       ref2    2       30      15=1X3=1X1=     *       0       0       GGTTTTATAAAACAAATAATT   ?????????????????????
x3      0       ref2    6       30      9=4I13= *       0       0       TTATAAAACAAATAATTAAGTCTACA      ??????????????????????????
x4      0       ref2    10      30      1X3=1X20=       *       0       0       CAAATAATTAAGTCTACAGAGCAAC       ?????????????????????????
x5      0       ref2    12      30      2=1X21= *       0       0       AATAATTAAGTCTACAGAGCAACT        ????????????????????????
x6      0       ref2    14      30      1X22=   *       0       0       TAATTAAGTCTACAGAGCAACTA ???????????????????????
```

### Usage in the literature

This tool was cited in Extensive sequencing of seven human genomes to characterize benchmark reference materials Sci Data. 2016; 3: 160025..

