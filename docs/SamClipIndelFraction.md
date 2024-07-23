# SamClipIndelFraction

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract clipping/indel fraction from BAM


## DEPRECATED

This tool can be replace with Bioalcidaejdk

## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samclipindelfraction  [options] Files

Usage: samclipindelfraction [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * clip



## Creation Date

20141118

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFractionTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFractionTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samclipindelfraction** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ samtools view -h -F3844 my.bam  | java -jar dist/jvarkit.jar samclipindelfraction

##UNMAPPED_READS=0
##MAPPED_READS=3028359
##CLIPPED_READS=1182730
##CLIPPED_READS_5_PRIME=597757
##CLIPPED_READS_3_PRIME=617399
##UNCLIPPED_READS=1845629
##COUNT_BASES=338644685
#CLIP	COUNT	FRACTION_OF_MAPPED_READS
0	1845629	0.5
1	7	1.8963724562195327E-6
2	6756	0.0018302703306027376
3	695	1.8828269386751074E-4
4	794	2.1510281860547272E-4
5	819	2.2187557737768533E-4
6	471	1.275987752684857E-4
7	447	1.210969268471616E-4
(...)
```

