# HmmMergeBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

For @154sns. Merge hmm bed


## Usage

```
Usage: hmmmergebed [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -n, --num
      dictionary if data are sorted on a REF order.
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference, --dict
      optional dictionary if data are sorted on a REF order.A SAM Sequence 
      dictionary source: it can be a *.dict file, a fasta file indexed with 
      'picard CreateSequenceDictionary', or any hts file containing a 
      dictionary (VCF, BAM, CRAM, intervals...)
    --version
      print version and exit
    -u
      Hide invalid lines
      Default: false

```


## Keywords

 * bed
 * merge


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew hmmmergebed
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20191108

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/hmm/HmmMergeBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/hmm/HmmMergeBed.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/hmm/HmmMergeBedTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/hmm/HmmMergeBedTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **hmmmergebed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a set of sorted BED file , or a file with the '.list' suffix containing the paths to the bed file. 

## Example:

```
$ head -n 3  ~/DC10* 
==> DC1074_UCSC.bed.txt <==
chr22	0	64000	6_Heterochrom_lowsignal
chr22	64000	65400	1_Polycomb_repressed
chr22	65400	119600	6_Heterochrom_lowsignal

==> DC1076_UCSC.bed.txt <==
chr22	0	119600	6_Heterochrom_lowsignal
chr22	119600	120000	3_Weak_promoter
chr22	120000	120200	4_Active_promoter

==> DC1077_UCSC.bed.txt <==
chr22	0	118400	7_Heterochrom_lowsignal
chr22	118400	119600	1_Polycomb_repressed
chr22	119600	120200	2_Poised_promoter

==> DC1082_UCSC.bed.txt <==
chr22	0	119600	6_Heterochrom_lowsignal
chr22	119600	120200	3_Weak_promoter
chr22	120200	122200	6_Heterochrom_lowsignal


$ java -jar dist/hmmmergebed.jar -n 3 ~/DC10*.txt 
chr22	0	64000	6_Heterochrom_lowsignal
#chr22	64000	65400	1_Polycomb_repressed:1;6_Heterochrom_lowsignal:2;7_Heterochrom_lowsignal:1
chr22	65400	119600	6_Heterochrom_lowsignal
chr22	119600	120000	3_Weak_promoter
#chr22	120000	120200	2_Poised_promoter:1;3_Weak_promoter:1;4_Active_promoter:2
#chr22	120200	120400	1_Polycomb_repressed:1;3_Weak_promoter:1;4_Active_promoter:1;6_Heterochrom_lowsignal:1
chr22	120400	122000	6_Heterochrom_lowsignal
chr22	122000	122600	3_Weak_promoter
#chr22	122600	122800	3_Weak_promoter:2;6_Heterochrom_lowsignal:2
chr22	122800	130000	6_Heterochrom_lowsignal
(...)

```

