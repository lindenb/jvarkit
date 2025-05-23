# BamStats04

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Coverage statistics for a BED file.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bamstats04  [options] Files

Usage: bamstats04 [options] Files
  Options:
  * -B, --bed
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
    -cov, --cov
      add this min coverage value to ask wether the position is not covered. 
      Use with care: any depth below this treshold will be trimmed to zero.
      Default: []
    -f, --filter, --jexl
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
    -partition, --partition
      [20171120]Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -R, --ref
      [20180126]If set, a column with the GC% will be added. Also used to read 
      CRAM. Indexed fasta Reference file. This file must be indexed with 
      samtools faidx and with picard/gatk CreateSequenceDictionary or samtools 
      dict 
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * coverage
 * depth
 * statistics
 * bed



## See also in Biostars

 * [https://www.biostars.org/p/309673](https://www.biostars.org/p/309673)
 * [https://www.biostars.org/p/348251](https://www.biostars.org/p/348251)



## Creation Date

20130513

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats04.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats04.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats04Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bamstats04/BamStats04Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamstats04** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## deprecated

use 'samtools coverage'

## input

Input is one or more indexed BAM files.

## History

* 2018-01-30: now using a jexl parser
* 2018-01-30: allow multiple values for '-cov'
* 2018-01-29: fixed bug from previous release (no data produced if no read). Added BioDas Resource.
* 2017-11-20: added new column 'partition'
* 2017-11-20: can read more than one BAM File.

## Example

```
$ java -jar dist/jvarkit.jar bamstats04 -B src/test/resources/toy.bed.gz src/test/resources/toy.bam 2> /dev/null | column -t 

#chrom  start  end  length  sample  mincov  maxcov  meancov  mediancov  nocoveragebp  percentcovered
ref     10     13   3       S1      3       3       3.0      3.0        0             100
ref2    1      2    1       S1      2       2       2.0      2.0        0             100
ref2    13     14   1       S1      6       6       6.0      6.0        0             100
ref2    16     17   1       S1      6       6       6.0      6.0        0             100
```

## Cited in:

  *  Han Ming Gan & al. , Genomic evidence of neo-sex chromosomes in the eastern yellow robin, GigaScience, Volume 8, Issue 9, September 2019, giz111, https://doi.org/10.1093/gigascience/giz111
  * Sigeman, H., Downing, P.A., Zhang, H. et al. The rate of W chromosome degeneration across multiple avian neo-sex chromosomes. Sci Rep 14, 16548 (2024). https://doi.org/10.1038/s41598-024-66470-7


