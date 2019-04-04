# BamCmpCoverage

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates the figure of a comparative view of the depths sample vs sample. Memory consideration: the tool alloc an array of bits which size is: (MIN(maxdepth-mindepth,pixel_width_for_one_sample) * count_samples)^2


## Usage

```
Usage: bamcmpcoverage [options] Files
  Options:
    -b, --bed
      restrict to region
    --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
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
    -M, --maxDepth
      max depth
      Default: 1000
    -m, --minDepth
      min depth
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    -r, --region
      restrict to region
    --version
      print version and exit
    -w, --width
      image width
      Default: 1000

```


## Keywords

 * sam
 * bam
 * visualization
 * coverage


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamcmpcoverage
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamCmpCoverage.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BamCmpCoverage.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamcmpcoverage** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)


```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```






### Screenshot

![img](https://pbs.twimg.com/media/B3in9wrIAAElLz8.jpg)

### Example

```
$ java -jar distBamCmpCoverage.jar  -o out.png file1.bam file2.bam fileN.bam
```


