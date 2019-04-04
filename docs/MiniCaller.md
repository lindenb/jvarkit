# MiniCaller

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Simple and Stupid Variant Caller designed for @AdrienLeger2


## Usage

```
Usage: minicaller [options] Files
  Options:
    -f, --filter
      [20171130](replaced with jexl expression). A JEXL Expression that will 
      be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -d, --mindepth
      Min depth
      Default: 20
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed Genome Reference. A fasta file that must be indexed with 
      samtools faidx and with picard CreateSequenceDictionary.
    -r, --region
      An interval as the following syntax : "chrom:start-end" or 
      "chrom:middle+extend"  or "chrom:start-end+extend" or 
      "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * calling
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew minicaller
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/calling/MiniCaller.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/calling/MiniCaller.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **minicaller** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Cited-In

  * "Direct Head-to-Head Evaluation of Recombinant Adeno-associated Viral Vectors Manufactured in Human versus Insect Cells". Kondratov & al. Molecular Therapy. [https://doi.org/10.1016/j.ymthe.2017.08.003](https://doi.org/10.1016/j.ymthe.2017.08.003).

## Example

```bash
$  java -jar dist/minicaller.jar -R ref.fa  bam.list > out.vcf
```

```
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Depth ReforAlt|Strand : RF,RR,AF,AR">
##FORMAT=<ID=DPG,Number=G,Type=Integer,Description="Depth for each allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is indel">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	4	.	T	A	.	.	DP=65	GT:DP:DP4:DPG	0/1:22:20,0,2,0:20,2	./.	./.	./.
rotavirus	5	.	T	A	.	.	DP=83	GT:DP:DP4:DPG	0:27:27,0,0,0:27	./.	0/1:25:20,0,5,0:20,5	./.
rotavirus	6	.	T	A	.	.	DP=97	GT:DP:DP4:DPG	0:30:30,0,0,0:30	0/1:21:20,0,1,0:20,1	0/1:33:31,0,2,0:31,2	./.
rotavirus	7	.	T	A	.	.	DP=112	GT:DP:DP4:DPG	0/1:38:36,0,2,0:36,2	0/1:23:21,0,2,0:21,2	0:37:37,0,0,0:37	./.
rotavirus	8	.	A	C	.	.	DP=122	GT:DP:DP4:DPG	0/1:41:38,0,3,0:38,3	0/1:26:25,0,1,0:25,1	0/1:40:38,0,2,0:38,2	./.
rotavirus	9	.	A	C	.	.	DP=139	GT:DP:DP4:DPG	0/1:46:44,0,2,0:44,2	0/1:29:27,0,2,0:27,2	0/1:48:44,0,4,0:44,4	./.
```

