# BamsToBedGraph

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

conver multiple BAM coverage to BedGraph


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bams2bedgraph  [options] Files

Usage: bams2bedgraph [options] Files
  Options:
    -B, --bed
      restrict to that input bed (or use --interval or reference dict)
    --depth
      if greater than 0, switch to mode : print the fraction of samples with a 
      depth greater or equal to 'x'
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --interval, --regions
      restrict to that interval chrom:start-end (or use --bed or reference 
      dict) 
    --mapq, -Q
      min mapping quality
      Default: 1
    --no-normalize
      do not divide by the number of bases (instead produces cumulative number 
      of bases)
      Default: false
    --output, -o
      Output file. Optional . Default: stdout
    --per_base
      output every base in each interval, do not group positions with the same 
      depth 
      Default: false
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * bam
 * alignment
 * bed
 * bedgraphs



## Creation Date

20260916

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2bedgraph/BamsToBedGraph.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2bedgraph/BamsToBedGraph.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2bedgraph/BamsToBedGraphTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2bedgraph/BamsToBedGraphTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bams2bedgraph** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


