# VcfDistanceBetweenVariants

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate variants with the distance between previous and next variant.


## Usage

```
Usage: vcfdistancevariants [options] Files
  Options:
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
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfdistancevariants
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190410

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfDistanceBetweenVariants.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfDistanceBetweenVariants.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfDistanceBetweenVariantsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfDistanceBetweenVariantsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfdistancevariants** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/vcfdistancevariants.jar src/test/resources/toy.vcf.gz
##fileformat=VCFv4.2
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=DIST_NEXT,Number=1,Type=Integer,Description="Distance to next variant">
##INFO=<ID=DIST_PREV,Number=1,Type=Integer,Description="Distance to previous variant">
(..)
##contig=<ID=ref,length=45>
##contig=<ID=ref2,length=40>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
ref	11	.	A	C	0.52	.	AC1=1;AF1=0.498013;BQB=1;DIST_NEXT=3;DP=3;DP4=2,0,1,0;FQ=-6.97257;MQ=30;MQ0F=0;MQB=1;PV4=1,1,1,0.106148;RPB=1;SGB=-0.379885	GT:PL	0/1:21,0,46
ref	14	.	AG	AGGG	0.71	.	AC1=2;AF1=1;DIST_PREV=3;DP=3;DP4=0,0,1,0;FQ=-37.5301;IDV=2;IMF=0.666667;INDEL;MQ=30;MQ0F=0;SGB=-0.379885	GT:PL	0/1:30,3,0
ref2	14	.	C	T	0.09	.	AC1=1;AF1=0.487148;BQB=1;DIST_NEXT=0;DP=6;DP4=3,0,1,0;FQ=-14.1557;MQ=30;MQ0F=0;MQB=1;PV4=1,0,1,0.0285955;RPB=1;SGB=-0.379885	GT:PL	0/0:13,0,64
ref2	14	.	CAA	CAAATAA	32.43	.	AC1=2;AF1=1;DIST_NEXT=1;DIST_PREV=0;DP=6;DP4=0,0,3,0;FQ=-43.5253;IDV=1;IMF=0.166667;INDEL;MQ=30;MQ0F=0;SGB=-0.511536;VDB=0.354794	GT:PL	1/1:72,9,0
ref2	17	.	T	A	0.14	.	AC1=1;AF1=0.491968;BQB=1;DIST_PREV=1;DP=6;DP4=4,0,1,0;FQ=-12.2521;MQ=30;MQ0F=0;MQB=1;PV4=1,1,1,0.201057;RPB=1;SGB=-0.379885	GT:PL	0/0:15,0,72
```

