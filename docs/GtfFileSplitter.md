# GtfFileSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split GTF file per gene, transcript, chromosome...


## Usage

```
Usage: gtfsplitter [options] Files
  Options:
    -compress, --compress, --gzip
      Gzip output gtf
      Default: false
    -H, --header
      include gtf header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -manifest, --manifest
      Manifest file containing the path to each gtf
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -m, --method
      split method. group: will create no more than 'pool-size' gtf files 
      where all the genes will be dispatched. stack: will create files where 
      there won't be more than 'pool-size' gene per file
      Default: gene
      Possible Values: [gene, transcript, contig, group, stack]
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    -p, --pool-size
      size for method=group or method=stack
      Default: 100
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * gtf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gtfsplitter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190807

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfFileSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfFileSplitter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf/GtfFileSplitterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf/GtfFileSplitterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtfsplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

split per gene

```
$ java -jar dist/gtfsplitter.jar -m gene -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"
TMP/51/4db9b1fbf7674369310d8698c987b6/22_ENSG00000100403.gtf
TMP/b1/315320e53be12e62a70acb6f26a2a6/3_ENSG00000183873.gtf
TMP/b1/035d795b2e064aa90cde9c64d17014/1_ENSG00000134250.gtf

$ cat TMP/b1/315320e53be12e62a70acb6f26a2a6/3_ENSG00000183873.gtf | cut -f 3 | sort | uniq -c 
    266 CDS
    282 exon
     17 five_prime_utr
      1 gene
     11 start_codon
     10 stop_codon
      9 three_prime_utr
     14 transcript

# check same number in the GTF file for this gene

$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep ENSG00000183873  | cut -f 3 | sort | uniq -c
    266 CDS
    282 exon
     17 five_prime_utr
      1 gene
     11 start_codon
     10 stop_codon
      9 three_prime_utr
     14 transcript
 
```

split per transcript:

```
$ java -jar dist/gtfsplitter.jar -m transcript -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"

TMP/7f/fd225f553e6fae72d87a26d21db3e6/3_ENST00000464652.gtf
TMP/ee/4b029d630d9d04ea90c5b6d58d6476/3_ENST00000455624.gtf
(...)
TMP/d9/e65eed61df7e80fdc9a12301b3b37d/22_ENST00000352645.gtf
TMP/c9/214537e6b2363128029fdd7898d1c7/1_ENST00000479412.gtf
TMP/23/750b00f3d4b5760ebaeaeda9e39e7f/3_ENST00000327956.gtf


$ cut -f 3 TMP/7f/fd225f553e6fae72d87a26d21db3e6/3_ENST00000464652.gtf
gene
exon
transcript
exon

# check the original file (the gene is ignored with this simple grep )
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep ENST00000464652 | cut -f 3 
exon
transcript
exon

```

split per chromosome/contig

```
$ java -jar dist/gtfsplitter.jar -m contig -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"
TMP/c4/ca4238a0b923820dcc509a6f75849b/1.gtf
TMP/b6/d767d2f8ed5d21a44b0e5886680cb9/22.gtf
TMP/ec/cbc87e4b5ce2fe28308fd9f2a7baf3/3.gtf

$ cut -f 1 TMP/b6/d767d2f8ed5d21a44b0e5886680cb9/22.gtf | uniq
22

```

per group

```
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep -v "#" | cut -f 3 | grep gene | wc
3

$ java -jar dist/gtfsplitter.jar -m group -p 2 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" -exec grep -w -c -H gene '{}' ';'
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:1
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf:2

# with a larger gtf file

$ java -jar dist/gtfsplitter.jar -m group -p 10 -o TMP ~/jeter.gtf.gz  
INFO	2019-08-07 11:14:00	SortingCollection	Creating merging iterator from 58 files

$ find TMP/ -name "*.gtf"    -exec grep -w -c gene -H '{}' ';'
TMP/ad/7732b2eb185c413c329c1f124bb809/Group000006.gtf:6229
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:6229
TMP/38/b153383e90aaa85a260d3560d36725/Group000005.gtf:6229
TMP/4b/6c5b9a6d1830f91089bdebe8321e90/Group000002.gtf:6230
TMP/2d/e7f75cf324bb4968ef397601e67fa7/Group000003.gtf:6229
TMP/93/2ea97d325f88f3c572d8af057b1c06/Group000008.gtf:6229
TMP/96/27eea1b15190b6a726b98af9f27b3a/Group000004.gtf:6229
TMP/96/ab60170db743f9704390147fadd705/Group000007.gtf:6229
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf:6230
TMP/45/a4780078f03085345ac98b5b2a766e/Group000009.gtf:6229

```

per stack

```
$ java -jar dist/gtfsplitter.jar -m stack -p 3 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" 
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
lindenb@mcclintock:~/src/jvarkit-git$ grep -w gene -c TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
3

$ java -jar dist/gtfsplitter.jar -m stack -p 1 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" 
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
TMP/4b/6c5b9a6d1830f91089bdebe8321e90/Group000002.gtf
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf

# with a larger gtf file

$ java -jar dist/gtfsplitter.jar -m stack -p 1000 -o TMP ~/jeter.gtf.gz 
INFO	2019-08-07 11:11:09	SortingCollection	Creating merging iterator from 58 files

$ find TMP/ -name "*.gtf"    -exec grep -w -c gene -H '{}' ';'
TMP/af/958f0b968e6c49f3474431cd3aa2df/Group000060.gtf:1000
TMP/ad/7732b2eb185c413c329c1f124bb809/Group000006.gtf:1000
TMP/e0/6a9eb20414e7406ac7a2309ebbe00f/Group000041.gtf:1000
TMP/e0/d2c4d002a926bea4295dbf6db3e831/Group000045.gtf:1000
TMP/f1/d98b194f295a688f1087400ef4cf6b/Group000020.gtf:1000
TMP/f1/07c7ca91843874df6a2c875ea674fd/Group000053.gtf:1000
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:1000
TMP/71/64a6dd2fa9469c1d475d9d6cb78cdc/Group000012.gtf:1000
TMP/a0/3524b44e90561115f6dd0ca3f5a00c/Group000016.gtf:1000
(...)

```



