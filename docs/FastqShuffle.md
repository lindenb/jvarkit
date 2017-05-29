# FastqShuffle

Shuffle Fastq files


## Usage

```
Usage: fastqshuffle [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -i
      single input is paired reads interleaved.(optional)
      Default: false
    -r
      random
      Default: java.util.Random@15316f2

```


## Keywords

 * fastq


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
$ make fastqshuffle
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqShuffle.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqShuffle.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqshuffle** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Synopsis

```
$ cat f.fq | java -jar dist/fastqshuffle.jar [options] 
$ java -jar dist/fastqshuffle.jar [options] f.fq.gz
$ java -jar dist/fastqshuffle.jar [options] f1.fq.gz f2.fq.gz

```


## Example

```bash
$ $ curl -s "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
java -jar dist/fastqshuffle.jar -i

@H06HDADXX130110:1:2103:11970:57672/1
GGATAGGGTTAGGGTTAGGGTTAGGGCTAGGGATAGGGGTAGGGTTGGGGTTGGTCATCGGGTGTTTCTTTGTGTTTGAGGTTGATTATTGTGATGGTTAAGGTATCTAGGTATTGTAAAAGTTGGCTTTTAACTTAGAAAATTATGTCATTCTGTTCACAAGTGTTTAGATTGGTAGATAGGTACTATGCGATCACTTCCATTGGCTGAGAGTTCGATTGATTATGAGCCACGCTAGTGGTTGAGATCT
+
69+26933-:7;;135,53<>7<692(?2=9:**;<=#####################################################################################################################################################################################################################
@H06HDADXX130110:1:2103:11970:57672/2
AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTATCGTCAAACCTTACCTCCTCCCTAGCCTCCACCCTGACCATGACACCAACCATCAGCCTTATAGAAAACCCCAGAGATGCTCTTATCCTATACCACAATTACCCCATAACGAAAGAAAGGACTGAAAACAAATAAGTAAAATTCGTACAAATTATATCTATGAGTATGTCCCTGAGTGTAGGTGTAGGTGCATCC
+
=>:=>@=?<>>??>;:<?<=;<<?>=;:8;=(5)0-6;1:>?<>##############################################################################################################################################################################################################
@H06JUADXX130110:1:1108:6424:55322/1
AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACTCTAACCCTAACCCTAACCCTAACGGTAACCCTTACCCTTACTGTAACGCTTATCCTAAATCAAATTCTTCCTCTTAAGATCGCTGTTAAAATTAATCCTATTAGAACAGGTCTTCTGGCACCAAGTTATGTCAATATCCCTTACTCTAAACATGCCTTGATCTCTCATGCATCACTTCAGCACAGCTCTTATGGATCTAGGATCCTCAGT
+
=>;=?=@@=?@?@@9>7@=?=;=?@>29?=?;=>@;4@*0878;40'=@;(3399@9>7@:A############################################################################################################################################################################################
@H06JUADXX130110:1:1108:6424:55322/2
AGGGATAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGATAGGGCTAGGGTTAGGGATAGGGATAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTATCGATAGGGATAGGGATAGGGATAGAGTTAGGGCTATGGGTAGGGTTAGAGTCAGGGAAAGAGATAGGGATGGAGATGGGGTTAAAAAGAAGTCAAGGAATTAAGGTAGGGAAACGGTTCGAGATCTGTAAAGGGCAACGA
+
>>;>*9?:@??@@????@????>@?>>@>@?>?????@@???????=<??8;*;:>?;+A?@?>89?@######################################################################################################################################################################################
@H06HDADXX130110:2:2116:3345:91806/1
GTTAGGGTTAGGGTTGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGCTAGGGTTAAGGGTAGGGTTAGCGAAAGGGCTGGGGTTAGGGGTGCGGGTACGCGTAGCATTAGGGCTAGAAGTAGGATCTGCAGTGCCTGACCGCGTCTGCGCGGCGACTGCCCAAAGCCTGGGGCCGACTCCAGGCTGAAGCTCAT
+
>=<=???>?>???=??>>8<?><=2=<===1194<?;:?>>?#3==>###########################################################################################################################################################################################################
@H06HDADXX130110:2:2116:3345:91806/2
TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCCTAACCCTAACCCTAACCCGTACCCTAAACCCAACCCTAACCACAAAGCAAATCCCAACCTTAACCGGAACCCGAAATCTCGCAGCAAATCTGCAGTAGAGACGCAGACTCAACCATGCGTCTATTAGTACGCATTATCATTGCCTCATGCTTCTTAAGTACAGAGAGATGAC
+
==;<?>@@@<>>@??<>>???<=>>?>:><@?4=:>7=5=>:<=@;'@A?########################################################################################################################################################################################################
```



