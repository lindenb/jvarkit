# Biostar170742

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert sam format to axt Format


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar biostar170742  [options] Files

Usage: biostar170742 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * axt



## See also in Biostars

 * [https://www.biostars.org/p/170742](https://www.biostars.org/p/170742)



## Creation Date

20151228

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar170742.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar170742.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar170742** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```
$ java -jar dist-2.0.1/biostar170742.jar \
	-R ref.fa \
	S1.bam  | head -n 15
 
0 rotavirus 1 70 rotavirus_1_317_5:0:0_7:0:0_2de/1 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT

1 rotavirus 1 70 rotavirus_1_535_4:0:0_4:0:0_1a6/2 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT

2 rotavirus 1 70 rotavirus_1_543_5:0:0_11:0:0_390/2 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT

3 rotavirus 1 70 rotavirus_1_578_3:0:0_7:0:0_7c/1 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT
(...)
```


