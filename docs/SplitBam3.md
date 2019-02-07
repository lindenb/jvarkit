# SplitBam3

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split a BAM by chromosome group


## Usage

```
Usage: splitbam3 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -g, --groupfile
      Chromosome group file. Interval are 1 based
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --mock
      add mock record if no samRecord saved in bam
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    -u, --unmapped
      unmapped chromosome name
      Default: Unmapped
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew splitbam3
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/splitbam/SplitBam3.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/splitbam/SplitBam3.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **splitbam3** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Split a BAM by chromosome group. Create EMPTY bams if no reads was found for a given group.
![img](https://chart.googleapis.com/chart?chl=+digraph+G+%7B%0D%0ABWA+-%3E+SPLITBAM%5Blabel%3D%22stdout%22%5D%3B%0D%0A+++SPLITBAM-%3ECHR1_bam%3B%0D%0A+++SPLITBAM-%3ECHR2_bam%3B%0D%0A+++SPLITBAM-%3ECHR3_bam%3B+%0D%0A+++CHR1_bam+-%3E+CHR1_vcf%3B%0D%0A+++CHR2_bam+-%3E+CHR2_vcf%3B%0D%0A+++CHR3_bam+-%3E+CHR3_vcf%3B%0D%0A+++CHR1_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR2_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR3_vcf+-%3E+merged_vcf%3B%0D%0A+%7D%0D%0A++++++++&cht=gv)

file output MUST contain the word '__GROUPID__'



### Example


the content of 'split_g1k_v37_01.txt'



```

CHROMS_01_09	1 2 3 4 5 6 7 8 9
CHROMS_10_0Y	10 11 12 13 14 15 16 17 18 19 20 21 22 X Y 
CHROMS_OTHER	MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 

```



split the output of bwa sampe on the fly:



```

bwa mem (...) | samtools sort (...) | \
java -jar dist/splitbam3.jar \
	-o TESTSPLITBAM/__GROUPID__.bam \
	-m \
	-g split_g1k_v37_01.txt 


[Fri Jul 26 13:25:56 CEST 2013] Executing as lindenb@master on Linux 2.6.32-358.6.2.el6.x86_64 amd64; OpenJDK 64-Bit Server VM 1.7.0_19-mockbuild_2013_04_17_19_18-b00; Picard version: null
INFO	2013-07-26 13:25:56	SplitBam	reading stdin
INFO	2013-07-26 13:25:56	SplitBam	opening TESTSPLITBAM/CHROMS_01_09.bam
INFO	2013-07-26 13:25:57	SplitBam	opening TESTSPLITBAM/CHROMS_10_0Y.bam
INFO	2013-07-26 13:25:58	SplitBam	opening TESTSPLITBAM/CHROMS_OTHER.bam
INFO	2013-07-26 13:35:58	SplitBam	closing group CHROMS_01_09
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_10_0Y
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_OTHER
INFO	2013-07-26 13:36:00	SplitBam	closing group Unmapped
Runtime.totalMemory()=1916600320

```







