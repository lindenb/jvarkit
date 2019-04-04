# VCFStripAnnotations

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Removes one or more field from the INFO/FORMAT column of a VCF.


## DEPRECATED

Use bcftools annotate -x 

## Usage

```
Usage: vcfstripannot [options] Files
  Options:
    -x, --exclude
      Use bcftools syntax INFO/x,INFO/y
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfstripannot
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstripannot/VCFStripAnnotations.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstripannot/VCFStripAnnotations.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstripannot** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ curl -sL "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
java -jar dist/vcfstripannot.jar -x "INFO/CSQ,INFO/EFF,INFO/AC,INFO/BaseQRankSum'  |\
grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M10475	M10478	M10500	M128215
chr1	145273345	.	T	C	289.85	.	AF=0.38;DP=1000;DS;Dels=0.00;FS=3.974;HRun=1;HaplotypeScore=17.4275;MQ=29.25;MQ0=0;MQRankSum=-1.370;QD=0.39;ReadPosRankSum=-1.117	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	156011444	.	T	C	2523.46	.	AF=0.50;DP=204;Dels=0.00;FS=4.328;HRun=0;HaplotypeScore=4.3777;MQ=35.24;MQ0=0;MQRankSum=-0.101;QD=14.93;ReadPosRankSum=1.575	GT:AD:DP:GQ:PL	0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	64982321	.	T	C	61.12	.	AF=1.00;DP=4;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=20.37	GT:AD:DP:GQ:PL	1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1142208	.	T	C	3404.30	.	AF=1.00;DP=122;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=2.6747;MQ=36.00;MQ0=0;QD=27.90	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126678092	.	G	A	89.08	.	AF=0.13;DP=185;Dels=0.00;FS=3.490;HRun=0;HaplotypeScore=3.3843;MQ=25.32;MQ0=0;MQRankSum=6.568;QD=2.02;ReadPosRankSum=-5.871	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135210791	.	T	C	65.41	.	AF=0.50;DP=11;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.2489;MQ=35.12;MQ0=0;MQRankSum=0.248;QD=16.35;ReadPosRankSum=-1.001	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	48873835	.	G	A	58.95	.	AF=1.00;DP=3;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=19.65	GT:AD:DP:GQ:PL	./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36779424	.	G	A	128.76	.	AF=0.13;DP=196;Dels=0.00;FS=1.447;HRun=0;HaplotypeScore=4.5749;MQ=36.22;MQ0=0;MQRankSum=-0.814;QD=3.90;ReadPosRankSum=-0.570	GT:AD:DP:GQ:PL	0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17819377	.	T	C	7515.25	.	AF=1.00;DP=319;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=7.7850;MQ=36.33;MQ0=0;QD=23.56	GT:AD:DP:GQ:PL	1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```

## History

* April 2017 : switched to BCFTOOLS syntax

