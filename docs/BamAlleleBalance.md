# BamAlleleBalance

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compute statistics about allele balance from a set of Bams


## Usage

```
Usage: bamallelebalance [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      Min MAPQ
      Default: 1
    --min-depth
      Ignore sites having a depth lower than 'x'
      Default: 10
    -o, --out
      Output file. Optional . Default: stdout
    --groupby, --partition
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -r, --range, --ratios
      Alleles Ratios. A 'range of double' is a list of floating number in 
      ascending order separated with semicolons.
      Default: com.github.lindenb.jvarkit.math.RangeOfDoubles@726f3b58
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
  * -v, --vcf, --variants
      Positions to consider. Must be frequence. Only unfiltered diallelic snp 
      are considered.
    --version
      print version and exit

```


## Keywords

 * vcf
 * allele-balance
 * depth
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamallelebalance
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200805

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/allelebalance/BamAlleleBalance.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/allelebalance/BamAlleleBalance.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamallelebalance** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Example

```bash
$ java -jar dist/bamallelebalance.jar --vcf src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S*.bam
RANGE	S5	S3	S2	S1	S4
[-Inf / 0.1[	0	1	1	1	1
[0.1 / 0.2[	0	0	0	0	0
[0.2 / 0.3[	0	0	0	0	0
[0.3 / 0.4[	0	0	0	0	0
[0.4 / 0.6[	0	0	0	1	1
[0.6 / 0.7[	0	0	0	0	1
[0.7 / 0.8[	0	0	0	0	2
[0.8 / 0.9[	0	0	0	0	0
[0.9 / Inf[	0	0	0	0	0
```

plotting the output

```R
T<-read.table("jeter.tsv", header=TRUE, sep="\t")
categories<-T[,1]
values<-as.matrix(T[,2:ncol(T)])
barplot(values, col = terrain.colors(length(categories)),legend=TRUE,xlab="Samples",ylab="Count",main="AD per Sample")
legend("topleft", 
       legend = categories, 
        bty = "n",
       fill = terrain.colors(length(categories))
    )
```

