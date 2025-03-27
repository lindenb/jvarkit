# FindAllCoverageAtPosition

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ?


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar findallcoverageatposition  [options] Files

Usage: findallcoverageatposition [options] Files
  Options:
    -clip, --clip
      use clipped bases (see also --extend).
      Default: false
    -x, --extend
      extend by 'x' base to try to catch close with clipped reads. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 500
    -filter, --filter
      [20171201](moved to jexl). A JEXL Expression that will be used to filter 
      out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-cigar
      Hide Cigar operators
      Default: false
    --left-align
      Left-aligns any indels in the read data contained in a BAM or CRAM file. 
      The same indel can often be placed at multiple positions and still 
      represent the same haplotype. While it is a commonly used convention to 
      place an indel at the left-most position
      Default: false
    -Q, --mapq
      Min mapping quality. Dicard reads having MAPQ < 'x'
      Default: 1
    -o, --out
      Output file. Optional . Default: stdout
    --groupby, --partition
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -f, --posfile
      File containing positions. if file suffix is '.bed': all positions in 
      the range will be scanned.
    -p, --position
      -p chrom:pos . Multiple separated by space. Add this chrom/position. 
      Required 
      Default: []
    -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict 
      Required for CRAM file or left-align
    --version
      print version and exit

```


## Keywords

 * bam
 * coverage
 * search
 * depth



## See also in Biostars

 * [https://www.biostars.org/p/259223](https://www.biostars.org/p/259223)
 * [https://www.biostars.org/p/250099](https://www.biostars.org/p/250099)
 * [https://www.biostars.org/p/409942](https://www.biostars.org/p/409942)
 * [https://www.biostars.org/p/9566474](https://www.biostars.org/p/9566474)



## Creation Date

20141128

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/findallcov/FindAllCoverageAtPosition.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/findallcov/FindAllCoverageAtPosition.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findallcoverageatposition** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Input

The input is a file containing a list of path to the bam.
 
## Example

```
$ find ./testdata/ -type f -name "*.bam" | \
 java -jar dist/findallcoverageatposition.jar -p rotavirus:100


#File              CHROM      POS  SAMPLE  DEPTH  M    I  D  N  S   H  P  EQ  X  Base(A)  Base(C)  Base(G)  Base(T)  Base(N)  Base(^)  Base(-)
./testdata/S4.bam  rotavirus  100  S4      126    126  0  0  0  29  0  0  0   0  5        0        0        121      0        0        0
./testdata/S1.bam  rotavirus  100  S1      317    317  1  0  0  50  0  0  0   0  27       0        1        289      0        1        0
./testdata/S2.bam  rotavirus  100  S2      311    311  0  1  0  60  0  0  0   0  29       1        0        281      0        0        1
./testdata/S3.bam  rotavirus  100  S3      446    446  1  0  0  86  0  0  0   0  39       0        1        406      0        1        0

```

## Cited in

 * Alexandra Bergfort, Tarek Hilal, Benno Kuropka, Ibrahim Avsar Ilik, Gert Weber, Tugce Aktas, Christian Freund, Markus C Wahl, The intrinsically disordered TSSC4 protein acts as a helicase inhibitor, placeholder and multi-interaction coordinator during snRNP assembly and recycling, Nucleic Acids Research, 2022;, gkac087, https://doi.org/10.1093/nar/gkac087
 * Tretter, C., de Andrade Kratzig, N., Pecoraro, M. et al. Proteogenomic analysis reveals RNA as a source for tumor-agnostic neoantigen identification. Nat Commun 14, 4632 (2023). https://doi.org/10.1038/s41467-023-39570-7

## See also

 * [https://twitter.com/pjacock/status/538300664334798848](https://twitter.com/pjacock/status/538300664334798848)
 * [https://twitter.com/yokofakun/status/538300434109456385](https://twitter.com/yokofakun/status/538300434109456385)
 * [https://twitter.com/pjacock/status/538299549455233024](https://twitter.com/pjacock/status/538299549455233024)
 * FindAVariation



