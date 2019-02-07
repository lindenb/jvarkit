# SamScanSplitReads

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

scan split reads


## Usage

```
Usage: samscansplitreads [options] Files
  Options:
    -B, --bed
      limit to that bed file
    --buffer-size
      dump buffer every 'x' bases. Most users should not use this.
      Default: 10000
    -x, --extend
      extends interval by 'x' pb before merging.
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -partition, --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -r, --rgn
      limit to that region CHROM:START-END
    --version
      print version and exit

```


## Keywords

 * sam
 * sv
 * splitreads
 * clip


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samscansplitreads
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReads.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReads.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReadsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReadsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samscansplitreads** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

finds the regions having some clipped reads.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
$ java -jar dist/samscansplitreads.jar src/test/resources/S*.bam 2> /dev/null 
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=M3,Number=1,Type=Float,Description="Median size of the clip in 3'">
##FORMAT=<ID=M5,Number=1,Type=Float,Description="Median size of the clip in 5'">
##FORMAT=<ID=N3,Number=1,Type=Integer,Description="Number of clipped reads in 3'">
##FORMAT=<ID=N5,Number=1,Type=Integer,Description="Number of clipped reads in 5'">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##contig=<ID=RF01,length=3302>
##contig=<ID=RF02,length=2687>
##contig=<ID=RF03,length=2592>
##contig=<ID=RF04,length=2362>
##contig=<ID=RF05,length=1579>
##contig=<ID=RF06,length=1356>
##contig=<ID=RF07,length=1074>
##contig=<ID=RF08,length=1059>
##contig=<ID=RF09,length=1062>
##contig=<ID=RF10,length=751>
##contig=<ID=RF11,length=666>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S3	S4	S5	S1	S2
RF01	195	.	N	<SPLIT>	.	.	DP=1;END=199;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:4.00:0:1	./.	./.
RF01	509	.	N	<SPLIT>	.	.	DP=1;END=577;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	725	.	N	<SPLIT>	.	.	DP=2;END=793;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF01	903	.	N	<SPLIT>	.	.	DP=2;END=1000;SVLEN=98	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	./.	./.
RF01	1607	.	N	<SPLIT>	.	.	DP=2;END=1616;SVLEN=10	GT:DP:M3:M5:N3:N5	0/1:1:0:9.00:0:1	./.	./.	./.	0/1:1:0:9.00:0:1
RF01	1672	.	N	<SPLIT>	.	.	DP=2;END=1682;SVLEN=11	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	./.	0/1:1:0:3.00:0:1	./.
RF01	1822	.	N	<SPLIT>	.	.	DP=1;END=1890;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	1926	.	N	<SPLIT>	.	.	DP=1;END=1994;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF01	2377	.	N	<SPLIT>	.	.	DP=2;END=2385;SVLEN=9	GT:DP:M3:M5:N3:N5	0/1:1:0:8.00:0:1	./.	./.	./.	0/1:1:0:8.00:0:1
RF01	2542	.	N	<SPLIT>	.	.	DP=2;END=2610;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:4.00:1:1	./.	./.
RF01	2689	.	N	<SPLIT>	.	.	DP=1;END=2691;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:2.00:0:1	./.	./.	./.
RF01	2719	.	N	<SPLIT>	.	.	DP=2;END=2787;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF01	3230	.	N	<SPLIT>	.	.	DP=2;END=3231;SVLEN=2	GT:DP:M3:M5:N3:N5	0/1:1:0:1.00:0:1	./.	./.	./.	0/1:1:0:1.00:0:1
RF02	3	.	N	<SPLIT>	.	.	DP=2;END=6;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF02	343	.	N	<SPLIT>	.	.	DP=4;END=451;SVLEN=109	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	0/1:2:68.00:0:2:0	0/1:1:0:3.00:0:1
RF02	513	.	N	<SPLIT>	.	.	DP=1;END=581;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	661	.	N	<SPLIT>	.	.	DP=1;END=729;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	818	.	N	<SPLIT>	.	.	DP=4;END=848;SVLEN=31	GT:DP:M3:M5:N3:N5	0/1:2:0:8.00:0:2	./.	./.	./.	0/1:2:0:8.00:0:2
RF02	957	.	N	<SPLIT>	.	.	DP=1;END=966;SVLEN=10	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	./.	./.	./.
RF02	1095	.	N	<SPLIT>	.	.	DP=1;END=1163;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF02	1707	.	N	<SPLIT>	.	.	DP=2;END=1725;SVLEN=19	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:9.00:0:1	0/1:1:0:2.00:0:1	./.	./.
RF02	1811	.	N	<SPLIT>	.	.	DP=1;END=1821;SVLEN=11	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:10.00:0:1	./.
RF02	1883	.	N	<SPLIT>	.	.	DP=2;END=1951;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF02	2220	.	N	<SPLIT>	.	.	DP=1;END=2224;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:4.00:0:1	./.
RF02	2515	.	N	<SPLIT>	.	.	DP=3;END=2663;SVLEN=149	GT:DP:M3:M5:N3:N5	./.	./.	0/1:3:68.00:4.00:2:1	./.	./.
RF03	500	.	N	<SPLIT>	.	.	DP=2;END=569;SVLEN=70	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	0/1:1:0:2.00:0:1	./.	./.
RF03	739	.	N	<SPLIT>	.	.	DP=2;END=823;SVLEN=85	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:2:68.00:2.00:1:1	./.
RF03	1072	.	N	<SPLIT>	.	.	DP=3;END=1147;SVLEN=76	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	0/1:1:68.00:0:1:0	./.	0/1:1:0:2.00:0:1
RF03	1207	.	N	<SPLIT>	.	.	DP=3;END=1350;SVLEN=144	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	0/1:1:68.00:0:1:0	./.
RF03	1729	.	N	<SPLIT>	.	.	DP=3;END=1809;SVLEN=81	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	0/1:1:0:7.00:0:1	./.	./.	0/1:1:68.00:0:1:0
RF03	1924	.	N	<SPLIT>	.	.	DP=2;END=1926;SVLEN=3	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	./.	./.	0/1:1:0:2.00:0:1
RF03	2153	.	N	<SPLIT>	.	.	DP=1;END=2221;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF04	173	.	N	<SPLIT>	.	.	DP=1;END=181;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:8.00:0:1	./.	./.	./.
RF04	579	.	N	<SPLIT>	.	.	DP=2;END=678;SVLEN=100	GT:DP:M3:M5:N3:N5	./.	./.	0/1:2:68.00:0:2:0	./.	./.
RF04	704	.	N	<SPLIT>	.	.	DP=2;END=707;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF04	754	.	N	<SPLIT>	.	.	DP=1;END=822;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF04	879	.	N	<SPLIT>	.	.	DP=1;END=887;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:8.00:0:1	./.
RF04	966	.	N	<SPLIT>	.	.	DP=2;END=1091;SVLEN=126	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	0/1:1:68.00:0:1:0	./.
RF04	1119	.	N	<SPLIT>	.	.	DP=1;END=1187;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF04	1378	.	N	<SPLIT>	.	.	DP=1;END=1380;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:2.00:0:1	./.	./.
RF04	1793	.	N	<SPLIT>	.	.	DP=2;END=1920;SVLEN=128	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:2:68.00:0:2:0	./.
RF04	2070	.	N	<SPLIT>	.	.	DP=1;END=2071;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:1.00:0:1	./.
RF05	112	.	N	<SPLIT>	.	.	DP=2;END=115;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF05	252	.	N	<SPLIT>	.	.	DP=1;END=256;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:0:4.00:0:1	./.
RF05	427	.	N	<SPLIT>	.	.	DP=1;END=433;SVLEN=7	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:6.00:0:1	./.	./.
RF05	529	.	N	<SPLIT>	.	.	DP=1;END=530;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:1.00:0:1	./.	./.
RF05	560	.	N	<SPLIT>	.	.	DP=2;END=563;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF05	750	.	N	<SPLIT>	.	.	DP=1;END=754;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	./.	./.
RF05	841	.	N	<SPLIT>	.	.	DP=1;END=909;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF05	960	.	N	<SPLIT>	.	.	DP=1;END=1028;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	./.	0/1:1:68.00:0:1:0	./.
RF05	1434	.	N	<SPLIT>	.	.	DP=1;END=1436;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:2.00:0:1	./.	./.	./.
RF06	26	.	N	<SPLIT>	.	.	DP=1;END=29;SVLEN=4	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:3.00:0:1	./.	./.
RF06	253	.	N	<SPLIT>	.	.	DP=1;END=321;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF06	465	.	N	<SPLIT>	.	.	DP=2;END=533;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF06	691	.	N	<SPLIT>	.	.	DP=2;END=762;SVLEN=72	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	0/1:1:68.00:0:1:0	./.
RF06	1045	.	N	<SPLIT>	.	.	DP=3;END=1133;SVLEN=89	GT:DP:M3:M5:N3:N5	./.	0/1:2:68.00:3.00:1:1	./.	0/1:1:68.00:0:1:0	./.
RF06	1224	.	N	<SPLIT>	.	.	DP=1;END=1292;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF07	223	.	N	<SPLIT>	.	.	DP=2;END=225;SVLEN=3	GT:DP:M3:M5:N3:N5	0/1:1:0:2.00:0:1	./.	./.	./.	0/1:1:0:2.00:0:1
RF07	345	.	N	<SPLIT>	.	.	DP=2;END=420;SVLEN=76	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	0/1:1:0:4.00:0:1	./.
RF07	790	.	N	<SPLIT>	.	.	DP=1;END=792;SVLEN=3	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:2.00:0:1	./.	./.
RF07	845	.	N	<SPLIT>	.	.	DP=2;END=913;SVLEN=69	GT:DP:M3:M5:N3:N5	0/1:1:68.00:0:1:0	./.	./.	./.	0/1:1:68.00:0:1:0
RF08	54	.	N	<SPLIT>	.	.	DP=2;END=57;SVLEN=4	GT:DP:M3:M5:N3:N5	0/1:1:0:3.00:0:1	./.	./.	./.	0/1:1:0:3.00:0:1
RF08	295	.	N	<SPLIT>	.	.	DP=1;END=296;SVLEN=2	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:0:1.00:0:1	./.	./.
RF08	668	.	N	<SPLIT>	.	.	DP=1;END=736;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	0/1:1:68.00:0:1:0	./.	./.	./.
RF08	896	.	N	<SPLIT>	.	.	DP=1;END=904;SVLEN=9	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:8.00:0:1	./.	./.	./.
RF10	192	.	N	<SPLIT>	.	.	DP=1;END=196;SVLEN=5	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:4.00:0:1	./.	./.	./.
RF10	433	.	N	<SPLIT>	.	.	DP=1;END=501;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
RF11	7	.	N	<SPLIT>	.	.	DP=3;END=133;SVLEN=127	GT:DP:M3:M5:N3:N5	./.	0/1:1:0:5.00:0:1	0/1:1:68.00:0:1:0	0/1:1:68.00:0:1:0	./.
RF11	179	.	N	<SPLIT>	.	.	DP=1;END=247;SVLEN=69	GT:DP:M3:M5:N3:N5	./.	./.	0/1:1:68.00:0:1:0	./.	./.
```
## history

 * 2018-11-14 rewritten from scratch

