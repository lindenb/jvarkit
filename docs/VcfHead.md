# VcfHead

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

print the first variants of a vcf


## Usage

```
Usage: vcfhead [options] Files
  Options:
    -c, --bycontig
      Print first variant for each contig; Implies VCF is sorted
      Default: false
    -n, -N, --count
      number of variants
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
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
$ ./gradlew vcfhead
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfHead.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfHead.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfHeadTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfHeadTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfhead** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 

```bash
$ curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
 java -jar dist/vcfhead.jar -n 2 | grep -v "##"

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  1094PC0005  1094PC0009  1094PC0012  1094PC0013
chr1    30860   .   G   C   33.46   .   AC=2;AF=0.053;AN=38;BaseQRankSum=2.327;DP=49;Dels=0.00;EFF=DOWNSTREAM(MODIFIER||||85|FAM138A|protein_coding|COD
ING|ENST00000417324|),DOWNSTREAM(MODIFIER|||||FAM138A|processed_transcript|CODING|ENST00000461467|),DOWNSTREAM(MODIFIER|||||MIR1302-10|miRNA|NON_CODING
|ENST00000408384|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000469289|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST000004
73358|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING
|ENST00000430492|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene
|NON_CODING|ENST00000488147|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000538476|);FS=3.128;HRun=0;HaplotypeScore=0.6718;In
breedingCoeff=0.1005;MQ=36.55;MQ0=0;MQRankSum=0.217;QD=16.73;ReadPosRankSum=2.017 GT:AD:DP:GQ:PL  0/0:7,0:7:15.04:0,15,177    0/0:2,0:2:3.01:0,3,39   0
/0:6,0:6:12.02:0,12,143    0/0:4,0:4:9.03:0,9,119
chr1    69270   .   A   G   2694.18 .   AC=40;AF=1.000;AN=40;DP=83;Dels=0.00;EFF=SYNONYMOUS_CODING(LOW|SILENT|tcA/tcG|S60|305|OR4F5|protein_coding|CODI
NG|ENST00000335137|exon_1_69091_70008);FS=0.000;HRun=0;HaplotypeScore=0.0000;InbreedingCoeff=-0.0598;MQ=31.06;MQ0=0;QD=32.86 GT:AD:DP:GQ:PL  ./. ./. 1/
1:0,3:3:9.03:106,9,0  1/1:0,6:6:18.05:203,18,0
```
 
