# BamXtremDepth

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compute low and high depth shared by a set of BAM files


## Usage

```
Usage: bamxtremdepth [options] Files
  Options:
    --contig
      limit to this chromosome.
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --interval_list
      force interval list output (default is BED)
      Default: false
    -Q, --mapq
      min mapping quality
      Default: 1
    --max-coverage, -M
      inclusive max coverage
      Default: 1000
    --min-coverage, -m
      inclusive min coverage
      Default: 0
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit
    -A
      alternative algorithm. Scan the whole chromosome instead of scanning the 
      previously found interval. Could be faster if there are too many 
      intervals. 
      Default: false
    -o
      Output file. Optional . Default: stdout

```


## Keywords

 * bam
 * depth
 * coverage


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamxtremdepth
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210511

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamxdepth/BamXtremDepth.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamxdepth/BamXtremDepth.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamxtremdepth** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of indexed BAM/CRAM files or a file with the '.list' suffix containing the paths to the BAM/CRAM paths.

## Example

```
$ java -jar dist/bamxtremdepth.jar --interval_list -M 15 -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null 
@HD	VN:1.6
@SQ	SN:RF01	LN:3302	M5:59dccb944425dd61f895a564ad7b56a7	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF02	LN:2687	M5:2c9c1ac1b7468ffaae96ad91c095c8b5	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF03	LN:2592	M5:d41f980f20d9cefbfd11ba2c1078f078	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF04	LN:2362	M5:935a2ad8d2f573d50b8c427f3b2f7c5d	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF05	LN:1579	M5:cdbaf0c352b44a79bef98daba1940d8a	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF06	LN:1356	M5:bdad6ae494bde13504d6988f2c7d94cd	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF07	LN:1074	M5:fa90eeb699fb6090b8f94611f19ac219	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF08	LN:1059	M5:4f8b3f714dc2327655941e6386c96b3b	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF09	LN:1062	M5:b221800f99aa2dc42258c0011ec228c8	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF10	LN:751	M5:f504f07ea3564b984207376aa02d8d00	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@SQ	SN:RF11	LN:666	M5:7a7cf2c7813f2e8bd74be383014202ca	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@CO	bamxtremdepth. compilation:20210511100604 githash:47c5dc0 htsjdk:2.23.0 date:20210511101013. cmd:--interval_list -M 15 -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
RF02	1	6	+	DP<=0
RF02	2686	2687	+	DP<=0
RF03	1	1	+	DP<=0
RF04	1	9	+	DP<=0
RF04	2358	2362	+	DP<=0
RF05	1	1	+	DP<=0
RF05	543	543	+	DP>=15
RF05	1573	1579	+	DP<=0
RF06	1	2	+	DP<=0
RF07	1	5	+	DP<=0
RF07	1074	1074	+	DP<=0
RF09	1	1	+	DP<=0
RF09	1062	1062	+	DP<=0
RF10	1	1	+	DP<=0
RF10	362	373	+	DP<=0
RF10	751	751	+	DP<=0
RF11	64	72	+	DP>=15
RF11	74	74	+	DP>=15
RF11	296	304	+	DP<=0
RF11	590	590	+	DP>=15
RF11	595	596	+	DP>=15

```


