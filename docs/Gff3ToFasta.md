# Gff3ToFasta

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

extract fasta from gtf


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gff2fasta  [options] Files

Usage: gff2fasta [options] Files
  Options:
    --fasta-length
      length of fasta lines
      Default: 60
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-cdna
      Hide cDNA
      Default: false
    --hide-mrna
      Hide mRNA
      Default: false
    --hide-peptide
      Hide peptide
      Default: false
    --hide-utr3
      Hide UTR3
      Default: false
    --hide-utr5
      Hide UTR5
      Default: false
    --region, --interval, --regions
      Limit to gene overlaping that region chr-end
    -o, --output
      Output file. Optional . Default: stdout
    --output-format
      output format
      Default: fasta
      Possible Values: [fasta, xml, bed]
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * gff
 * gff3
 * fasta



## Creation Date

20241016

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gff2fa/Gff3ToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gff2fa/Gff3ToFasta.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gff2fasta** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


convert GFF3 to fasta/xml/bed

## Example
```
$ java -jar dist/jvarkit.jar gff2fasta -R ref.fasta in.gff3 --output-format xml | xmllint --format -
(...)
  <entry>
    <attributes>
      <attribute key="gene_name">SCN5A</attribute>
      <attribute key="gff_type">transcript</attribute>
      <attribute key="start">38674799</attribute>
      <attribute key="length">138</attribute>
      <attribute key="type">UTR5</attribute>
      <attribute key="chrom">chr3</attribute>
      <attribute key="transcript_id">ENST00000327956.6</attribute>
      <attribute key="strand">negative</attribute>
      <attribute key="protein_id">ENSP00000333674.6</attribute>
      <attribute key="build">GRCh37</attribute>
      <attribute key="end">38687267</attribute>
      <attribute key="location">chr3:38674799-38687267</attribute>
      <attribute key="transcript_type">protein_coding</attribute>
      <attribute key="gene_id">ENSG00000183873.11</attribute>
    </attributes>
    <sequence>AGTGGACACTGTGGGCATGCGTGTCCGTGCAGCACATCGCCATGCAGGAGCCTGGGGGAAGGCCTTTCCCTCAGTGAGGGCTGCAGCTTCCCCACAGGCAACGTGAGGAGAGCCTGTGCCCAGAAGCAGGATGAGAAG</sequence>
  </entry>
</gff2fasta>

 java -jar dist/jvarkit.jar gff2fasta -R  ref.fasta in.gff3 --output-format bed | head | cut -c1-100 | column -t
chrom  start     end       strand    location                gene_id             gene_name  gff_type    length  transcript_id  transcript_name  tran
chr3   38589548  38674839  negative  chr3:38589548-38674840  ENSG00000183873.11  SCN5A      transcript  8303    ENST                            
chr3   38591812  38674797  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  5997    ENST                            
chr3   38591812  38674797  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  1999    ENST                            
chr3   38589548  38591810  negative  chr3:38589548-38591811  ENSG00000183873.11  SCN5A      transcript  2264    ENST                            
chr3   38674799  38674839  negative  chr3:38674799-38674840  ENSG00000183873.11  SCN5A      transcript  42      ENST00                          
chr3   38589553  38674852  negative  chr3:38589553-38674853  ENSG00000183873.11  SCN5A      transcript  8362    ENST                            
chr3   38591812  38674797  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  6048    ENST                            
chr3   38591812  38674797  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  2016    ENST                            
chr3   38589553  38591810  negative  chr3:38589553-38591811  ENSG00000183873.11  SCN5A      transcript  2259    ENST                            



```

