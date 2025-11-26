# KnownGeneToFasta

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert ucsc genpred to fasta


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar kg2fa  [options] Files

Usage: kg2fa [options] Files
  Options:
    --coding
      ignore non-coding transcripts.
      Default: false
    --empty
      Discard empty sequences.
      Default: false
    --hide, --exclude
      Exclude the following type of sequence: mRNA, cDNA, peptide, utr5, utr3 
      , uORF, uPeptide (case insensitive, comma/space separated)
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    -sql, --sql
      SQL Schema URI. Each instance of transcript can be associated to a .sql 
      schema to help the software to decode the semantics of the columns. Eg.: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV20.sql
      Default: <empty string>
    --version
      print version and exit
    -L
      fasta line length.
      Default: 50

```


## Keywords

 * kg
 * knownGene
 * fasta
 * genpred



## Creation Date

20190213

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2fa/KnownGeneToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2fa/KnownGeneToFasta.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2fa** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar kg2fa -R human_g1k_v37.fasta -D  --case 3 -L 0 | cut -c 1-200 | head

>ENST00000456328.2 585|ENST00000456328.2|chr1|+|11868|14409|11868|11868|3|11868,12612,13220,|12227,12721,14409,|0|DDX11L1|none|none|-1,-1,-1,
GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGG
>ENST00000607096.1 585|ENST00000607096.1|chr1|+|30365|30503|30365|30365|1|30365,|30503,|0|MIR1302-11|none|none|-1,
GGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTA
>ENST00000417324.1 585|ENST00000417324.1|chr1|-|34553|36081|34553|34553|3|34553,35276,35720,|35174,35481,36081,|0|FAM138A|none|none|-1,-1,-1,
CACACAACGGGGTTTCGGGGCTGTGGACCCTGTGCCAGGAAAGGAAGGGCGCAGCTCCTGCAATGCGGAGCAGCCAGGGCAGTGGGCACCAGGCTTTAGCCTCCCTTTCTCACCCTACAGAGGGCAGGCCCTTCAGCTCCATTCTCCTCCAAGGCTGCAGAGGGGGCAGGAATTGGGGGTGACAGGAGAGCTGTAAGGTC
>ENST00000335137.3 585|ENST00000335137.3|chr1|+|69090|70008|69090|70008|1|69090,|70008,|0|OR4F5|cmpl|cmpl|0,
ATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTCCTATTTATGTTGTTTTTTGTATTCTATGGAGGAATCGTGTTTGGAAACCTTCTTATTGTCATAACAGTGGTATCTGACTCCCACCTTCACTCTCCCATGTACTTCCTGCTAGCCAACCTCTCACTCATTGATCTGTCTCTGTC
>ENST00000466430.1 585|ENST00000466430.1|chr1|-|89294|120932|89294|89294|4|89294,92090,112699,120774,|91629,92240,112804,120932,|0|RP11-34P13.7|none|none|-1,-1,-1,-1,
CTGATCCATATGAATTCCTCTTATTAAGAAAAATAAAGCATCCAGGATTCAATGAAGAACTGACTATCACCTTGTTAATCATTCAGAAACATGTTGCAGGCTTAAGCCATTTTTGATATAGATACTGAAACAATTACTTGCTAAGAGCAAACTTGAAGgtatggataaggccctgagtcatcttcctgagctgaatgata
(...)

```
 

