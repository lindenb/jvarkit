# SamStructVarNeg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find Structural Variation by Negative Comparaison


## Usage

```
Usage: svneg [options] Files
  Options:
    -B, --bed
      Optional BED file. SV should overlap this bed.
    --chrom-regex
      Only consider the chromosomes matching the following regular expression.
      Default: (chr)?([1-9][0-9]*|X|Y)
  * -b, --controls, --bams
      Control Bams. One path per lines
    -d, --distance
      Max distance between two read to test if they both end at the same ~ 
      position. A distance specified as a positive integer.Comma are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --low-complexity
      zone of low complexity (=many discordant reads). In those region the 
      software will 'freeze'. If the number of buffered reads * this number 
      per sample then clear the buffer
      Default: 2000
    --mapq
      min mapping quality.
      Default: 0
    --max-contigs-in-buffer
      clear buffer if it contains more than 'x' different mate contig to 
      prevent area matching everywhere
      Default: 3
    --max-controls
      Maximum number controls bam matching the event
      Default: 1
    --max-sa
      ignore reads having more that SA:Z supplementary alignments.
      Default: 4
    -m, --min
      Min number of events to validate the translocation
      Default: 3
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -D
      Presence of a discordant read in the control, whatever is the contig or 
      the distance to theoritical mate is enought to invalidate the candiate
      Default: false

```


## Keywords

 * sam
 * bam
 * sv
 * translocation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew svneg
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190413

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNeg.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNeg.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNegTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamStructVarNegTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **svneg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



```
$ java -jar dist/svneg.jar --bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam

##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SR,Number=4,Type=Integer,Description="Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="other chromosome">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="other position">
##INFO=<ID=STDDEV_POS1,Number=1,Type=Integer,Description="std deviation to position 1">
##INFO=<ID=STDDEV_POS2,Number=1,Type=Integer,Description="std deviation to position 2">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variation type">
##svneg.meta=compilation:20190415102101 githash:4c6799f7 htsjdk:2.19.0 date:20190415102110 cmd:--bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG02260
9	137230996	9:137230996:14:79839048	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=14;DP=18;POS2=79839048;STDDEV_POS1=120;STDDEV_POS2=187;SVTYPE=BND	GT:DP:SR	0/1:18:5,13,13,5
14	79839119	14:79839119:9:137230968	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=9;DP=18;POS2=137230968;STDDEV_POS1=154;STDDEV_POS2=145;SVTYPE=BND	GT:DP:SR	0/1:18:13,5,5,13


```

