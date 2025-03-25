# SamViewWithMate

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract reads within given region(s), and their mates


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samviewwithmate  [options] Files

Usage: samviewwithmate [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
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
  * -b, --bed, -r, --region
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
      Default: (unspecified)
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    -st, --streaming
      Force Streaming mode even if bam is indexed. Warning: Streaming mode 
      doesn't garantee that all mates will be fetched because a read only 
      contains the start position of the mate of which may be out of the 
      user's intervals, unless the MC (mate cigar) attribute is defined.
      Default: false
    --supplementary
      Look also in the supplementary data of the SA:Z attribute
      Default: false
    -u, --unmapped
      Also search for the unmapped mates. Not available in streaming mode.
      Default: false
    --version
      print version and exit

```


## Keywords

 * sam
 * bam



## See also in Biostars

 * [https://www.biostars.org/p/151403](https://www.biostars.org/p/151403)
 * [https://www.biostars.org/p/105714](https://www.biostars.org/p/105714)
 * [https://www.biostars.org/p/368754](https://www.biostars.org/p/368754)
 * [https://www.biostars.org/p/9537698](https://www.biostars.org/p/9537698)



## Creation Date

20190207

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/viewmate/SamViewWithMate.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/viewmate/SamViewWithMate.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/viewmate/SamViewWithMateTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/viewmate/SamViewWithMateTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samviewwithmate** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Update 2025

you should also have a look at `samtools view --fetch-pair`

## How it works

Two modes:

  * The streaming mode is set if the input is `stdin` or if the bam file is NOT indexed. The input is scanned and any read (or mate) that overlap a region is written.

  * The other mode use the bam index. First we scan the regions, we collect the other regions and the names of the reads (requires memory), the bam is opened a second time and we collect the reads.

## Example

```
$ java -jar dist/samviewwithmate.jar -r "9:137230721-137230796"  ./src/test/resources/HG02260.transloc.chr9.14.bam | cut -f 1-9 | tail
ERR251239.10989793	83	9	137230747	60	30S70M	=	137230326	-490
ERR251239.3385449	147	9	137230754	60	1S99M	=	137230352	-500
ERR251240.17111373	99	9	137230764	60	100M	=	137231150	475
ERR251240.46859433	147	9	137230777	60	65S35M	=	137230342	-469
ERR251240.74563730	147	9	137230787	60	1S99M	=	137230407	-478
ERR251240.1291708	83	9	137230789	60	100M	=	137230411	-477
ERR251240.11887757	97	9	137230795	37	100M	14	79839451	0
ERR251239.34016218	81	14	79839349	37	100M	9	137230679	0
ERR251240.10196873	81	14	79839368	37	100M	9	137230721	0
ERR251240.11887757	145	14	79839451	37	100M	9	137230795	0
```

## Cited in

  * Garsed, D.W., Pandey, A., Fereday, S. et al. The genomic and immune landscape of long-term survivors of high-grade serous ovarian cancer. Nat Genet (2022). https://doi.org/10.1038/s41588-022-01230-9


