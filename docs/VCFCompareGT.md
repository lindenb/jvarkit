# VCFCompareGT

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

 compare two or more genotype-callers for the same individuals. Produce a VCF with FORMAT fields indicating if a genotype is new or modified.


## Usage

```
Usage: vcfcomparegt [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -label, --labels
      A comma separated list of label that will be used as the title of the 
      vcfs. Must be provided in the same order. If blank, some numeric indexes 
      will be used
      Default: <empty string>
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -nc, --nocall2homref
      convert no call to hom-ref
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -vf, --variant-filter
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    --version
      print version and exit
    -m
      only print modified samples
      Default: false

```


## Keywords

 * vcf
 * compare


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcomparegt
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGT.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGT.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGTTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFCompareGTTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparegt** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of VCF files or a file with '.list' suffix with the path (one path per line).

Genotypes are supposed diploids.

## Example

```bash

$ java -jar dist/vcfcomparegt.jar -m  Sample.samtools.vcf.gz Sample.gatk.vcf.gz

##fileformat=VCFv4.1
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GCH,Number=1,Type=Integer,Description="Changed Genotype">
##FORMAT=<ID=GNW,Number=1,Type=Integer,Description="Genotype Created/Deleted">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Qual">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=GDF,Number=.,Type=String,Description="Samples with Genotype Difference">
##VCFCompareGT_1=File: Sample.samtools.vcf.gz
##VCFCompareGT_2=File: Sample.gatk.vcf.gz
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_2	Sample_1
X	1860854	rs5781	A	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1866893	rs2824	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1878904	.	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	0/1:20:0:1:71	./.
X	1895117	.	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/0:2:0:1:27
X	1895755	.	C	AG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:4:0:1:17
X	1900009	rs6181	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:13:0:1:30	./.
X	1905130	.	AG	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:3:0:1:16
X	1905160	.	A	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1905165	.	C	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
X	1913889	.	C	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1948846	rs6	T	TG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:239:0:1:99	./.
X	1955199	.	C	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
(...)
```

