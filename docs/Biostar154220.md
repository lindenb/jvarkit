# Biostar154220

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Cap BAM to a given coverage


## Usage

```
Usage: biostar154220 [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -d, -n, --depth
      expected coverage.
      Default: 20
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
    --keep-unmapped
      write unmapped reads
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --query-sorted
      Input was sorted on query name but I promess there is one and only one 
      chromosome: e.g: samtools view -h in.bam 'chr1:234-567' | samtools sort 
      -n -) .
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * coverage
 * depth



## See also in Biostars

 * [https://www.biostars.org/p/154220](https://www.biostars.org/p/154220)
 * [https://www.biostars.org/p/9471266](https://www.biostars.org/p/9471266)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar154220
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20150812

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar154220Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar154220** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$ java -jar dist/sortsamrefname.jar --samoutputformat BAM input.bam  |\
  java -jar dist/biostar154220.jar  -n 20 --samoutputformat BAM |\
  samtools sort -T tmp -o output.bam -


$ samtools mpileup output.bam  | cut -f 4 | sort | uniq -c

  12692 0
 596893 1
  94956 10
  56715 11
  76947 12
  57912 13
  66585 14
  51961 15
  63184 16
  47360 17
  65189 18
  65014 19
 364524 2
 169064 20
  72078 3
 118288 4
  54802 5
  82555 6
  53175 7
  78474 8
  54052 9

```

## Cited in

  * "Burden of Cardiomyopathic Genetic Variation in Lethal Pediatric Myocarditis". Amy R. Kontorovich & al. 6 Jul 2021 https://doi.org/10.1161/CIRCGEN.121.003426Circulation: Genomic and Precision Medicine. 2021;0


