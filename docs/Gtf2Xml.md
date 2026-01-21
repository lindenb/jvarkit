# Gtf2Xml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert GTF/GFF to XML


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gtf2xml  [options] Files

Usage: gtf2xml [options] Files
  Options:
    -d, --distance
      if >=0, add a 'y' attribute that could be used to display the bed 
      records in a browser, 	this 'y' is the graphical row where the item 
      should be displayed. 	This distance is the distance between two item 
      where there is a collision. Memory consuming . A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --interval, --regions, --region, -R
      only in the following interval CHR:START-END. First Scan for genes in 
      the interval and extend the interval to get the whole genes. Requires 
      GTF input is a tabix index file.
      Default: <empty string>
    -o, --output
      Output file. Optional . Default: stdout
    --simple
      Don't print group data by gene/transcript. Print each GTF record on the 
      fly 
      Default: false
    --skip-attributes
      Don't print the following attributes. Multiple separated by 
      commas/spaces/semicolons 
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * xml
 * gtf
 * gff
 * gff3



## See also in Biostars

 * [https://www.biostars.org/p/478242](https://www.biostars.org/p/478242)



## Creation Date

20150811

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf2xml/Gtf2Xml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf2xml/Gtf2Xml.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf2xml/Gtf2XmlTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf2xml/Gtf2XmlTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtf2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$ java -jar dist/jvarkit gtf2xml src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | xmllint --format - | head -n 100
<?xml version="1.0" encoding="UTF-8"?>
<gtf genebuild-last-updated="2013-09" genome-build="GRCh37.p13" genome-build-accession="NCBI:GCA_000001405.14" genome-date="2009-02" genome-version="GRCh37">
  <gene id="ENSG00000100403" chrom="22" start="41697526" end="41756151" strand="+" source="ensembl_havana" type="gene">
    <attributes>
      <attribute key="gene_id">ENSG00000100403</attribute>
      <attribute key="gene_version">10</attribute>
      <attribute key="gene_name">ZC3H7B</attribute>
      <attribute key="gene_source">ensembl_havana</attribute>
      <attribute key="gene_biotype">protein_coding</attribute>
    </attributes>
    <transcripts>
      <transcript id="ENST00000486331" chrom="22" start="41697719" end="41732847" strand="+" source="havana" type="transcript">
        <attributes>
          <attribute key="gene_id">ENSG00000100403</attribute>
          <attribute key="gene_version">10</attribute>
          <attribute key="transcript_id">ENST00000486331</attribute>
          <attribute key="transcript_version">1</attribute>
          <attribute key="gene_name">ZC3H7B</attribute>
          <attribute key="gene_source">ensembl_havana</attribute>
          <attribute key="gene_biotype">protein_coding</attribute>
          <attribute key="transcript_name">ZC3H7B-002</attribute>
          <attribute key="transcript_source">havana</attribute>
          <attribute key="transcript_biotype">retained_intron</attribute>
          <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
          <attribute key="havana_transcript_version">1</attribute>
        </attributes>
        <exon chrom="22" start="41697719" end="41697776" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">1</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00001942555</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        <transcript chrom="22" start="41697719" end="41732847" strand="+" source="havana" type="transcript">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
          </attributes>
        </transcript>
        <exon chrom="22" start="41716659" end="41716717" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">2</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00003530265</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        <exon chrom="22" start="41721568" end="41721601" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">3</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00003553644</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        (...)
```


