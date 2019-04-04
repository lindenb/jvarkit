# VCFAnnoBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with the Coverage statistics of a BAM file+  BED file of capture. 


## DEPRECATED

useless: use DP/DP4 in the Genotypes, or use GATK variant annotator

## Usage

```
Usage: vcfannobam [options] Files
  Options:
    -BAM, --bam
      Indexed BAM File. One file ending with the '.list' suffix will be 
      interpreted as a text file containing the path to the bams.
  * -BED, --bed
      BED File Exome capture.
    -MIN_COV, --coverage
      min coverage to say the position is not covered
      Default: 0
    -filter, --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -tag, --tag
      VCF info tag
      Default: CAPTURE
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * depth
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfannobam
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannobam/VCFAnnoBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannobam/VCFAnnoBam.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfannobam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$  java -jar dist/vcfannobam.jar \
		-BAM input.bam\
		-BED capture.bed \
		input.vcf.gz

(...)
##INFO=<ID=CAPTURE,Number=1,Type=String,Description="Capture stats: Format is (start|end|mean|min|max|length|not_covered|percent_covered) ">
(...)
2	16100665	.	A	T	13261.77	.	CAPTURE=16100619|16100715|1331.96|1026.0|1773.0|97|0|100
2	178395141	.	T	A	1940.77	.	CAPTURE=178394991|178395199|193.11|100.0|276.0|209|0|100
(...)
```

