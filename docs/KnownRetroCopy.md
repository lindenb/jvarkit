# KnownRetroCopy

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate VCF structural variants that could be intron from retrocopies.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar knownretrocopy  [options] Files

Usage: knownretrocopy [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --distance
      max distance between an intron and the deletion found in the VCF
      Default: 10
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -k, --known
      Gene-ID of known retrogenes. One per line. A source could be : 
      http://retrogenedb.amu.edu.pl/static/download/ 
    -o, --out
      Output file. Optional . Default: stdout
  * -genpred, --genpred, --kg, --transcripts
      Transcrips as genpred format 
      https://genome.ucsc.edu/FAQ/FAQformat.html#format9  . The genePred 
      format is a compact alternative to GFF/GTF because one transcript is 
      described using only one line.	Beware chromosome names are formatted the 
      same as your REFERENCE. A typical KnownGene file is http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV47.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
    --version
      print version and exit

```


## Keywords

 * gtf
 * retrocopy
 * deletion



## Creation Date

20190815

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/KnownRetroCopy.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/KnownRetroCopy.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/KnownRetroCopyTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/KnownRetroCopyTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **knownretrocopy** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/knownretrocopy.jar --gtf Homo_sapiens.GRCh37.87.gtf.gz candidateSV.vcf.gz | grep RETR

##FILTER=<ID=RETROCOPY_INTRON,Description="variant could be a deleted intron from a retrocopy">
##FILTER=<ID=RETROCOPY_KNOWN,Description="variant could be a deleted intron from a known retrocopy">
##INFO=<ID=RETROCOPY,Number=.,Type=String,Description="Identifiers for the retrocopies.">
chr1	38077349	MantaDEL:2204:0:0:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
chr1	38077349	MantaDEL:2204:0:1:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
chr1	38077349	MantaDEL:2204:1:1:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
(...)
```


