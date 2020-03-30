# CoverageServer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Jetty Based http server serving Bam coverage.


## Usage

```
Usage: coverageserver [options] Files
  Options:
    --bed, --B
      Optional bed file containing user's intervals. 4th column is used as the 
      name of the interval
    -o, --output, --comment
      Output file for writing comments as a BED file. Very basic= not suitable 
      for multiple users.
    --extend
      Extend interval by this factor. e.g: if x='0.5' chr1:100-200 -> 
      chr1:50-250 
      Default: 1.0
    --gtf
      Optional Tabix indexed GTF file. Will be used to retrieve an interval by 
      gene name, or to display gene names in a region.
    --height
      Image height
      Default: 300
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --link, --url, --hyperlink
      creates a hyperlink when 'click' in an area. The URL must contains 
      __CHROM__, __START__ and __END__ that will be replaced by their values. 
      Predefined values are 'hg19','hg38','igv'. IGV : 
      "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" , 
      UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: <empty string>
    --images-per-row, --ipr
      Number of images per row.
      Default: 2
    --vcf, --region, --regions, --intervals
      Same as --bed but intervals won't be annotated. A source of intervals. 
      The following suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, 
      gff, gff.gz, gtf.gz.Otherwise it could be an empty string (no interval) 
      or a list of plain interval separated by '[ \t\n;,]'
      Default: (empty)
    --known
      Optional Tabix indexed Bed or VCF file containing known CNV. Both types 
      must be indexed.
    --mapq
      Min. Read Mapping Quality.
      Default: 0
    --max_-size
      Security for memory. Max interval size. A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb
      Default: 10000000
    --pedigree, -p
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --port
      server port.
      Default: 8080
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --sashimi
      Enable sashimi plots.
      Default: false
    --small-length
      show reads if the region has a length <= 'x'.
      Default: 1000
    --version
      print version and exit
    --width
      Image width
      Default: 500

```


## Keywords

 * cnv
 * bam
 * coverage
 * server


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew coverageserver
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200212

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/server/CoverageServer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/server/CoverageServer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **coverageserver** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a set of indexed BAM file or a file with the suffix `.list` containing the path to the bams.
 
## Example

```
java -jar dist/coverageserver.jar \
	--pedigree fam.ped \
	--bed roi.bed \
	-o comments.bed \
	-R ref.fasta src/test/resources/S*.bam

```
## Hidden http parameters

 * `columns=5` change the number of columns at runtime.

## Screenshot

![https://twitter.com/yokofakun/status/1227932501747871745](https://pbs.twimg.com/media/EQp-Ga4XsAAxNYn?format=png&name=small)

![https://twitter.com/yokofakun/status/1228260742157209601](https://pbs.twimg.com/media/EQuooeGX0AAAHeu?format=jpg&name=medium)

![https://twitter.com/yokofakun/status/1229343426036076546](https://pbs.twimg.com/media/EQ-BJSXWkAItBtJ?format=jpg&name=medium)

![https://twitter.com/yokofakun/status/1238112128646733824](https://pbs.twimg.com/media/ES6oQbmWoAAxAx9?format=png&name=small)


