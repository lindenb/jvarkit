# ShiftBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

shit all coordinates of a bam


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar shiftbam  [options] Files

Usage: shiftbam [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
  * -R2, --destination-reference
      Original fasta reference. We shift the bam back to this reference.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    -R, --reference, --source-reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --stringency
      Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * bed



## Creation Date

20241001

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/shiftbam/ShiftBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/shiftbam/ShiftBam.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **shiftbam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


shift coordinates of a BAM.

Say the bam has been mapped on a sub-fasta.

```
# here all contigs look like "chrxxx:123-456"
samtools faidx ref.fa "chr1:2345-6789" >  ref2.fa
samtools faidx ref.fa "chr21:12345-16789" >>  ref2.fa




(...)
bwa mem ref2.fa read1.fq read2.fq |\
	java -jar dist/jvarkit.jar shiftbam -R ref.Fa

````


In this sub-fasta all the chromosomes NAMES **MUST** look like 'chr1:123-456'.
This tool takes as input the original REF and the bam and it's coordinates in the new ref
and shit the coordinates to the new reference, including the SA:Z attributes.
Ouput BAM is not sorted on coordinate.



