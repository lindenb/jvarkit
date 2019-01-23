# GcPercentAndDepth

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extracts GC% and depth for multiple bam using a sliding window


## Usage

```
Usage: gcpercentanddepth [options] Files
  Options:
    -filter, --filter
      [20171219]A JEXL Expression that will be used to filter out some 
      sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -partition, --partition
      [20171219]Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    --version
      print version and exit
    -B
       (file.bed) (optional). If not defined: use whole genome. Warning memory 
      consumming: must alloc sizeof(int)*win.size()*num(samples).
    -R
      Indexed Genome Reference. It can be a the path to fasta file that must 
      be indexed with samtools faidx and with picard CreateSequenceDictionary. 
      It can also be a BioDAS dsn url like 
      `http://genome.cse.ucsc.edu/cgi-bin/das/hg19/` . BioDAS references are 
      slower, but allow to work without a local reference file.
    -m
       min depth
      Default: 0
    -n
       skip window if Reference contains one 'N'.
      Default: false
    -o
      Output file. Optional . Default: stdout
    -s
       (window shift)
      Default: 50
    -w
       (window size)
      Default: 100
    -x
       don't print genomic index.
      Default: false

```


## Keywords

 * gc%
 * depth
 * coverage


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gcpercentanddepth
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/GcPercentAndDepth.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/GcPercentAndDepth.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gcpercentanddepth** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## History

* 2017-12-09 : currently rewriting everything... do not use...

## Example

```bash
$ java -jar dist/gcanddepth.jar -R ref.fasta -b capture.bed 1.bam 2.bam ... > result.tsv
```
