# CommBams

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Equivalent of unix 'comm' for bams sorted on queryname


## Usage

```
Usage: commbams [options] Files
  Options:
    -delim, --delimiter
      Output delimiter
      Default: 	
    -empty, --empty
      Empty content symbol
      Default: .
    -f, --format
      What should I print ? (only the read name ? etc...)
      Default: name
      Possible Values: [name, but_metadata, all]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -1, --hide1
      suppress read unique to file 1
      Default: false
    -2, --hide2
      suppress read unique to file 2
      Default: false
    -3, --hide3
      suppress reads present in both files
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -sortmethod, --sortmethod
      [20171110]Method used to sort the read on query name. (samtools != 
      picard) see https://github.com/samtools/hts-specs/issues/5
      Default: picard
      Possible Values: [samtools, picard]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * comm
 * compare


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew commbams
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBams.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBams.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBamsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cmpbams/CommBamsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **commbams** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

### Example 1 

``` 
 $ java -jar dist/commbams.jar --samtools \
 	 -f but_metadata -delim '\n' \
 	B00GWFP_std.hg19.qname.bam B00GWFP_S1.hg19.qname.bam
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1	83	chr5	21564864	40	151M	=	21564540	-475	CTCCCAGAGAGAAGCATCAACAGCTTAGGGTGTAGTCTAAACAGAAATCTTGCACTCCTCCTGCAGTAGCGTCTCTATTTTTTATGCTGAACATTATTTGCTAATTCCAACTGGCTCTAAGCTAATGTGTTTCCCAGGTTTTCTCAATGAN	AFAA<,,,<,,,,,,,,7,,,,7,,,,A7KKF<,F,,7,7,A,A7F7,K<A,,,,7,,,7KKKFAFA,7,A7F7,7,,,KFF,,AKKFFFF<<K<<KAFAKA,A,,A7,AAAFKFA,A,FKKAA,AKKKKFFFKF<KKKKKKKFFAAAA<#
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1	77	*	0	0	*	*	0	0	NTCATTGAGAAAACCTGGGAAACACATTAGCTTAGAGCCAGTTGGAATTAGCAAATAATGTTCAGCATAAAAAATAGAGACGCTACTGCAGGAGGAGTGCAAGATTTCTGTTTAGACTACACCCTAAGCTGTTGATGCTTCTCTCTGGGA	!<AAAAFFKKKKKKK<FKFFFKKKKA,AAKKF,A,AFKFAAA,7A,,A,AKAFAK<<K<<FFFFKKA,,FFK,,,7,7F7A,7,AFAFKKK7,,,7,,,,A<K,7F7A,A,7,7,,F,<FKK7A,,,,7,,,,7,,,,,,,,<,,,<AAF
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2	163	chr5	21564540	60	8S106M37S	=	21564864	475	NTAAGAATATTTCACACTTAAAACAAAATCTGATTAGACAAACACTTTGATTGTTATTATTCGCGTATATCATCTACCAGAAGCAAATAGACATCTACTACATCTTTCAAGAAAGTTTACCTATCAATATTACTCAACTGGACCCAATAAT	#<A,<,,A,,K<7FKFF,7,,AF,,7,7AAF,,7<<,,7,AF,,7,7A<,7FA,,7,,7F,,A7FKK7,7,,,,,7,,,,,<,<,,,,7,<,,,,,,7FF7AF<7,,<,,,,7,7,,,,,,,<,,,,,,,,,,,,,,,,,,<,,,,,,,,,
PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2	141	*	0	0	*	*	0	0	NTAAGAATATTTCACACTTAAAACAAAATCTGATTAGACAAACACTTTGATTGTTATTATTCGCGTATATCATCTACCAGAAGCAAATAGACATCTACTACATCTTTCAAGAAAGTTTACCTATCAATATTACTCAACTGGACCCAATAA	!<A,<,,A,,K<7FKFF,7,,AF,,7,7AAF,,7<<,,7,AF,,7,7A<,7FA,,7,,7F,,A7FKK7,7,,,,,7,,,,,<,<,,,,7,<,,,,,,7FF7AF<7,,<,,,,7,7,,,,,,,<,,,,,,,,,,,,,,,,,,<,,,,,,,,
```

### Example 2 

```
 $ java -jar dist/commbams.jar --samtools \
 	B00GWFP_std.hg19.qname.bam B00GWFP_S1.hg19.qname.bam

.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:5388/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:6513/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:6513/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:7181/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:7181/2
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10205/1
.	.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10205/2
.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10380/1	.
.	PANORAMIX:1:HJY2CCCXX:7:1101:1133:10380/2	.
```


