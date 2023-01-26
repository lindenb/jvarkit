# GtfToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert GTF/GFF3 to BED.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gtf2bed  [options] Files

Usage: gtf2bed [options] Files
  Options:
    -c, --columns
      comma separated columns to be displayed
      Default: gtf.source,gtf.feature,gtf.score,gtf.strand,gtf.frame,ccds_id,exon_id,exon_number,exon_version,gene_biotype,gene_id,gene_name,gene_source,gene_version,havana_transcript,havana_transcript_version,protein_id,protein_version,tag,transcript_biotype,transcript_id,transcript_name,transcript_source,transcript_version,Alias,biotype,ccdsid,constitutive,description,ensembl_end_phase,ensembl_phase,exon_id,external_name,gene_id,havana_transcript,havana_version,ID,logic_name,Name,Parent,protein_id,rank,tag,transcript_id,version
    --grep
      Check some identifiers are found in a column. syntax: 
      <COLUMN>:<FILE_CONTAINING_THE_IDENTIFIERS> 
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * gtf
 * gff
 * gff3
 * bed



## Creation Date

20220629

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtf2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```
$ java -jar ${JVARKIT_DIST}/gtf2bed.jar --columns gtf.feature,gene_name,gene_biotype,gene_id ~/src/jvarkit/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | head | column -t
#chrom  start      end        gtf.feature      gene_name  gene_biotype    gene_id
1       120454175  120459317  exon             NOTCH2     protein_coding  ENSG00000134250
1       120454175  120612240  gene             NOTCH2     protein_coding  ENSG00000134250
1       120454175  120457928  three_prime_utr  NOTCH2     protein_coding  ENSG00000134250
1       120454175  120612240  transcript       NOTCH2     protein_coding  ENSG00000134250
1       120457928  120457931  stop_codon       NOTCH2     protein_coding  ENSG00000134250
1       120457931  120459317  CDS              NOTCH2     protein_coding  ENSG00000134250
1       120460287  120460385  CDS              NOTCH2     protein_coding  ENSG00000134250
1       120460287  120460385  exon             NOTCH2     protein_coding  ENSG00000134250
1       120461028  120461176  CDS              NOTCH2     protein_coding  ENSG00000134250
```


