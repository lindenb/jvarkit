# KnownDeletion

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find Split reads in a region to validate a known CNV


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar knowndel  [options] Files

Usage: knowndel [options] Files
  Options:
    -B, --bam
      Optional Bam to sasve the matching reads
    -x, --extend
      search in the boundaries of the CNV around +/- 'x' bases. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 100
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax. 'default' is 'mapqlt(1) || Duplicate() || 
      FailsVendorQuality() || NotPrimaryAlignment() || 
      SupplementaryAlignment()' 
      Default: mapqlt(1) || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      For Reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
  * -r, --rgn, --region
      region chr:start-end to be validated.
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * sam
 * vcf
 * sv
 * deletion



## Creation Date

20190523

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/KnownDeletion.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/KnownDeletion.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **knowndel** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

search known deletion
## Example

```
find DIR -type f -name "*.bam" > bam.list

```


