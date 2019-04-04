# VcfClusteredReadEdge

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Variant annotation :   variants clustered near the ends of reads


## Usage

```
Usage: vcfclusteredreadedge [options] Files
  Options:
    -B, --bams
      path of indexed BAM path with read Groups. You can put those paths in a 
      text file having a *.list sufffix
      Default: []
    -d, --distance
      minimal distance to the end of the **CLIPPED** read.
      Default: 1
    -filter, --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -gt, --gt
      Genotype FILTER name
      Default: EDGEVAR
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -vt, --vt
      Variant FILTER name: set if ALL Genotypes have a variant near the edge.
      Default: EDGEVAR

```


## Keywords

 * sam
 * bam
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfclusteredreadedge
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfClusteredReadEdge.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfClusteredReadEdge.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfClusteredReadEdgeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfClusteredReadEdgeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfclusteredreadedge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Pour Sandro B.

GATK ClusteredReadPosition https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_ClusteredReadPosition.php  only works with Mutect2

The program looks for SNV in the VCF, go back to the reads in the bam.

For one variant , if all the reads contain the variant at less than 'distance' then the genotype is FILTERED

if all the reads are FILTERED, the variant is FILTERED


## Example

```
java -jar dist/vcfclusteredreadedge.jar -B in.bam in.vcf
```

```
find . -name "*.bam" > tmp.list 
java -jar dist/vcfclusteredreadedge.jar -B tmp.list  in.vcf
```


