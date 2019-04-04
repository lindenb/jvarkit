# SamExtractClip

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract Soft Clipped Sequences from a SAM. Ouput is a FASTQ


## Usage

```
Usage: samextractclip [options] Files
  Options:
    -c, --clipped
      Print the original Read where the clipped regions have been removed.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --minsize
      Min size of clipped read
      Default: 5
    -p, --original
      Print Original whole Read that contained a clipped region.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -readFilter, --readFilter
      [20181208]A JEXL Expression that will be used to filter out some 
      sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * fastq
 * clip



## See also in Biostars

 * [https://www.biostars.org/p/125874](https://www.biostars.org/p/125874)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samextractclip
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamExtractClip.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamExtractClip.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamExtractClipTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamExtractClipTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samextractclip** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$ curl -L -s "https://raw.githubusercontent.com/samtools/samtools/develop/test/dat/mpileup.1.sam" |\
	java -jar dist/samextractclip.jar 2> /dev/null 

@ERR013140.3521432/1:0:17:1:99
AGAGGTCCCCAACTTCTTTGCA
+
@AEDGBHIIIIIFJGIKHGHIJ
@ERR156632.12704932/2:0:17:1:163
TGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
BFAFGFEIGFEFHHEIDKJGHHHJIIE=@KKGGKJG
@ERR156632.9601178/1:0:17:1:99
CTATGACAGGGAGGTCATGTGCAGGCTGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
DEEEIIHHKIJILKHLHIKEKHHMKLKKJGKKKKLKLFIHEKIKL=KLJLKIILHKMH9LJJ
@ERR162872.21706338/1:0:17:1:99
CTTCTTTGCA
+
BHBFH<EIFG
@ERR243039.1049231/2:0:17:1:163
TGCAGGCTGGAGAAGGGGACAAGAGGTCCCCAACTTCTTTGCA
+
AEEFIFHIJDGIGIJHHIAGGGLGJIEJHJHHFIJGJJDFJIG
@ERR013140.20277196/2:1:17:97:163
AAACTGTCCAGCGCATACCCGCATCCCCCCAAAAGAAGCCACCGCCCCAACACACACCCCCCACCCGCATAACC
+
00($,+3(*+..,%%+6%*#%2,/001)%%$2%%/$.%$00(,%+,1'*.%7(%&$&#'$$$#%#%#($+%+)"
@ERR013140.19887184/1:0:17:99:113
GTGTGTGTCGGGGGTGTCTGGGGTCTCACCCACGACCAAC
+
%$($&$*+#%%#1'$$%2-'0&3$/$/$-73/69:7=1%2
@ERR013140.4461044/1:0:17:114:113
ACTCCCTGGGCCTGGCA
+
/=1:/=44-348<0(91
@ERR013140.3521432/2:0:17:226:147
CACCCCTAGAAGTGACGGC
+
71%??A9A792/7-2%(&:
@ERR013140.11659627/1:0:17:645:83
TTAGCAACAAAAAGGAC
+
%5?-$)89<=;9>(.14
@ERR013140.7259970/1:0:17:660:83
ACGCCTGGTACA
+
40&/81&8:/<<
@ERR013140.29762488/1:0:17:716:83
GGACTCA
+
)/4/142
@ERR013140.11567710/1:0:17:984:83
TGCTTGA
+
/36>+5/
```


### History

* 20190221 : handle sam records without quality https://github.com/lindenb/jvarkit/issues/121
* 20180412 : fastq is now reverse complemented if read was on negative strand


### Cited In

 * Perlman syndrome nuclease DIS3L2 controls cytoplasmic non-coding RNAs and provides surveillance pathway for maturing snRNAs : http://nar.oxfordjournals.org/content/44/21/10437.full


