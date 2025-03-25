# KgToGff

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert UCSC genpred/knowngene file to gff3 or gtf


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar kg2gff  [options] Files

Usage: kg2gff [options] Files
  Options:
    --coding
      select coding transcript only.
      Default: false
    --gtf
      use the GTF writer instead of GFF3
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --modulo3
      discard transcripts where CDS length isn't a modulo 3 (eg. remove 
      xeno-transcript mapped on another build)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --score
      default score (-1: undefined)
      Default: -1
    --source
      label for column 'source' in the output
      Default: ucsc
    --version
      print version and exit

```


## Keywords

 * gff
 * gff3
 * ucsc
 * genpred



## See also in Biostars

 * [https://www.biostars.org/p/9610527](https://www.biostars.org/p/9610527)



## Creation Date

20210106

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2gff/KgToGff.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2gff/KgToGff.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2gff** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Warning

ouput is escaped for UTF8, some characters like ':' might be converted to a hexadecimal encoding.

# Example

```
headjeter.kg -n1 | java -jar dist/jvarkit.jar kg2gff 
##gff-version 3.1.25
chr1	ucsc	gene	55039456	55064852	.	+	.	ID=gene%3APCSK9;Name=PCSK9;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;gene_type=protein_coding
chr1	ucsc	mRNA	55039456	55064852	.	+	.	ID=ENST00000710286.1;Parent=gene%3APCSK9;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	exon	55039456	55040044	.	+	.	ID=ENST00000710286.1%3AE0;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE0
chr1	ucsc	exon	55043843	55044034	.	+	.	ID=ENST00000710286.1%3AE1;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE1
chr1	ucsc	exon	55046523	55046646	.	+	.	ID=ENST00000710286.1%3AE2;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE2
chr1	ucsc	exon	55052278	55052411	.	+	.	ID=ENST00000710286.1%3AE3;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE3
chr1	ucsc	exon	55052650	55052791	.	+	.	ID=ENST00000710286.1%3AE4;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE4
chr1	ucsc	exon	55055993	55056189	.	+	.	ID=ENST00000710286.1%3AE5;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE5
chr1	ucsc	exon	55057331	55057514	.	+	.	ID=ENST00000710286.1%3AE6;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE6
chr1	ucsc	exon	55058036	55058209	.	+	.	ID=ENST00000710286.1%3AE7;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE7
chr1	ucsc	exon	55058499	55058647	.	+	.	ID=ENST00000710286.1%3AE8;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE8
chr1	ucsc	exon	55059486	55059663	.	+	.	ID=ENST00000710286.1%3AE9;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE9
chr1	ucsc	exon	55061375	55061556	.	+	.	ID=ENST00000710286.1%3AE10;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE10
chr1	ucsc	exon	55063369	55064852	.	+	.	ID=ENST00000710286.1%3AE11;Parent=ENST00000710286.1;Name=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1;exon_id=ENST00000710286.1%3AE11
chr1	ucsc	CDS	55039481	55040044	.	+	0	ID=ENST00000710286.1%3AC1;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55043843	55044034	.	+	0	ID=ENST00000710286.1%3AC2;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55046523	55046646	.	+	0	ID=ENST00000710286.1%3AC3;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55052278	55052411	.	+	2	ID=ENST00000710286.1%3AC4;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55052650	55052791	.	+	0	ID=ENST00000710286.1%3AC5;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55055993	55056189	.	+	2	ID=ENST00000710286.1%3AC6;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55057331	55057514	.	+	0	ID=ENST00000710286.1%3AC7;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55058036	55058209	.	+	2	ID=ENST00000710286.1%3AC8;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55058499	55058647	.	+	2	ID=ENST00000710286.1%3AC9;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55059486	55059663	.	+	0	ID=ENST00000710286.1%3AC10;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55061375	55061556	.	+	2	ID=ENST00000710286.1%3AC11;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	CDS	55063369	55063584	.	+	0	ID=ENST00000710286.1%3AC12;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	5_prime_UTR	55039456	55039480	.	+	.	ID=ENST00000710286.1%3AU55039456;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
chr1	ucsc	3_prime_UTR	55063585	55064852	.	+	.	ID=ENST00000710286.1%3AU55063585;Parent=ENST00000710286.1;biotype=protein_coding;gene_id=gene%3APCSK9;gene_name=PCSK9;transcript_id=ENST00000710286.1
```


