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


