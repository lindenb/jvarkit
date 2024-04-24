# LowResBam2Raster

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Low Resolution BAM to raster graphics


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar lowresbam2raster  [options] Files

Usage: lowresbam2raster [options] Files
  Options:
    -clip, --clip
      Show clipping
      Default: false
    -depth, --depth
      Depth track height.
      Default: 100
    -gcPercent, --gcPercent
      GC% track height.
      Default: 100
    -gcwin, --gcWindowSize
      GC% Window size
      Default: 10
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -gtf, --gtf
      A GTF (General Transfer Format) file. See 
      https://www.ensembl.org/info/website/upload/gff.html . Please note that 
      CDS are only detected if a start and stop codons are defined.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hideInsert, --hideInsertions
      Hide insertions
      Default: false
    --highlight
      hightligth those positions.
      Default: []
    --mapqopacity
      How to handle the MAPQ/ opacity of the reads. all_opaque: no opacity, 
      handler 1: transparency under MAPQ=60
      Default: handler1
      Possible Values: [all_opaque, handler1]
    --limit, --maxrows
      Limit number of rows to 'N' lines. negative: no limit.
      Default: -1
    -minh, --minh
      Min. distance between two reads.
      Default: 10
    -noSuppl, --noSuppl
      Hide arcs of Supplementary alignments.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout [20180829] filename can be also 
      an existing directory or a zip file, in witch case, each individual will 
      be saved in the zip/dir.
    -printNames, --printNames
      Print Read Names (for debugging)
      Default: false
    -proper, --proper
      Hide read if in a paired-end pair, both reads are mapped but not in 
      proper pair.
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
  * -r, --region
      Restrict to that region. An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or "chrom:start-end+extend" 
      or "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
    -srf, --samRecordFilter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax. 'default' is 'mapqlt(1) || Duplicate() || 
      FailsVendorQuality() || NotPrimaryAlignment() || 
      SupplementaryAlignment()' 
      Default: mapqlt(1) || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    --single
      Convert paired reads to single-end reads.
      Default: false
    --spaceyfeature
      number of pixels between features
      Default: 1
    -V, --variants, --vcf
      VCF files used to fill the position to hightlight with POS
      Default: []
    --version
      print version and exit
    -w, --width
      Image width
      Default: 1000
    -D
      set some css style elements. '-Dkey=value'. Undocumented.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * png
 * gtf



## See also in Biostars

 * [https://www.biostars.org/p/293741](https://www.biostars.org/p/293741)



## Creation Date

20170523

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2Raster.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2Raster.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2RasterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2RasterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **lowresbam2raster** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/lowresbam2raster.jar \
	-o out.png -r "22:38999+10000" in.bam \
	 -clip -srf "" -R ref.fasta  -kg knownGene.txt.gz
```

## Cited in

 * Baudic M, Murata H, Bosada FM, Melo US, Aizawa T, Lindenbaum P, van der Maarel LE, Guedon A, Baron E, Fremy E, Foucal A, Ishikawa T, Ushinohama H, Jurgens SJ, Choi SH, Kyndt F, Le Scouarnec S, Wakker V, Thollet A, Rajalu A, Takaki T, Ohno S, Shimizu W, Horie M, Kimura T, Ellinor PT, Petit F, Dulac Y, Bru P, Boland A, Deleuze JF, Redon R, Le Marec H, Le Tourneau T, Gourraud JB, Yoshida Y, Makita N, Vieyres C, Makiyama T, Mundlos S, Christoffels VM, Probst V, Schott JJ, Barc J. TAD boundary deletion causes PITX2-related cardiac electrical and structural defects. Nat Commun. 2024 Apr 20;15(1):3380. doi: 10.1038/s41467-024-47739-x. PMID: 38643172; PMCID: PMC11032321.

## see also

  * https://twitter.com/yokofakun/status/951769190884610051
  * https://twitter.com/yokofakun/status/973836167522279425
  * https://twitter.com/notSoJunkDNA/status/1012309599079272448

## Screenshots

![https://pbs.twimg.com/media/DAldDxvXkAAGMoJ.jpg](https://pbs.twimg.com/media/DAldDxvXkAAGMoJ.jpg)

![https://pbs.twimg.com/media/DTVcmGYW4AAmiZp.jpg](https://pbs.twimg.com/media/DTVcmGYW4AAmiZp.jpg)

![https://pbs.twimg.com/media/DYPC0XFWAAAk5SV.jpg](https://pbs.twimg.com/media/DYPC0XFWAAAk5SV.jpg)

![https://pbs.twimg.com/media/Dgxp_5OXkAEbAYW.jpg](https://pbs.twimg.com/media/Dgxp_5OXkAEbAYW.jpg)


