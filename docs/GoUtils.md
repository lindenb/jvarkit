# GoUtils

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Gene Ontology Utils. Retrieves terms from Gene Ontology


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar goutils  [options] Files

Usage: goutils [options] Files
  Options:
    -A, --accession
      User Go Terms accession numbers or name.eg GO:0005216 ('ion channel 
      activity') 
      Default: []
    -af, --accession-file
      File containing accession numbers. One per line. After the first white 
      space one can define optional attributes for 
      gexf:`color=<COLOR>;size=<SIZE> 
    -action, --action
      What shoud I do ? default is dump as table. 'goa' only keeps GOA 
      elements in GOA input in GAF format (e.g 
      http://geneontology.org/gene-associations/goa_human.gaf.gz). 
      Default: dump_table
      Possible Values: [dump_table, dump_gexf, goa, gff3]
    --exclude-accession, --exclude
      User Go Terms to be EXCLUDED accession numbers or name.eg
      Default: []
    -g, --gff, --gff3
      GFF3 file for action=goa or action=gff3.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --inverse
      inverse the result
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -go
      Gene ontology URI. Formatted as OBO format.
      Default: http://current.geneontology.org/ontology/go-basic.obo
    -goa
      Gene Ontology Annotation GOA input in GAF format (e.g 
      http://geneontology.org/gene-associations/goa_human.gaf.gz). 
      Default: http://geneontology.org/gene-associations/goa_human.gaf.gz

```


## Keywords

 * geneontology
 * go
 * gexf



## See also in Biostars

 * [https://www.biostars.org/p/488538](https://www.biostars.org/p/488538)



## Creation Date

20180130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/go/GoUtils.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/go/GoUtils.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **goutils** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

children of  GO:0005216 'ion channel activity' 

```
$ java -jar dist/goutils.jar -R is_a  -A 'GO:0005216' 

#ACN	NAME	DEFINITION
GO:1905030	voltage-gated ion channel activity involved in regulation of postsynaptic membrane potential	Any voltage-gated ion channel activity that is involved in regulation of postsynaptic membrane potential.
GO:1905057	voltage-gated calcium channel activity involved in regulation of postsynaptic cytosolic calcium levels	Any voltage-gated calcium channel activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:1905054	calcium-induced calcium release activity involved in regulation of presynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of presynaptic cytosolic calcium ion concentration.
GO:1905058	calcium-induced calcium release activity involved in regulation of postsynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:0016286	small conductance calcium-activated potassium channel activity	Enables the transmembrane transfer of potassium by a channel with a unit conductance of 2 to 20 picoSiemens that opens in response to stimulus by internal calcium ions. Small conductance calcium-activated potassium channels are more sensitive to calcium than are large conductance calcium-activated potassium channels. Transport by a channel involves catalysis of facilitated diffusion of a solute (by an energy-independent process) involving passage through a transmembrane aqueous pore or channel, without evidence for a carrier-mediated mechanism.
GO:0043855	cyclic nucleotide-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0043854	cyclic nucleotide-gated mechanosensitive ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens in response to a mechanical stress and when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0099142	intracellularly ATP-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when ATP has been bound by the channel complex or one of its constituent parts on the intracellular side of the plasma membrane.
GO:0099101	G-protein gated potassium channel activity	A potassium channel activity that is gated by binding of a G-protein beta-gamma dimer.
(...)
```

## Example

### action = goa

```
$ wget -q -O - "http://geneontology.org/gene-associations/goa_human.gaf.gz" | gunzip -c | java -jar dist/goutils.jar -go go.obo --action goa -A 'ion transmembrane transport'  | head
UniProtKB	A0A1W2PN81	CHRNA7	enables	GO:0022848	GO_REF:0000002	IEA	InterPro:IPR002394	FNeuronal acetylcholine receptor subunit alpha-7	CHRNA7	protein	taxon:9606	20210612	InterPro
UniProtKB	A0PJK1	SLC5A10	enables	GO:0015370	Reactome:R-HSA-8876283	TAS		F	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20200515	Reactome
UniProtKB	A0PJK1	SLC5A10	involved_in	GO:0035725	GO_REF:0000108	IEA	GO:0015370	P	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903401	GO_REF:0000108	IEA	GO:0015189	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903826	GO_REF:0000108	IEA	GO:0015181	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	PMID:22946059	IMP		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20140424	UniProt
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	Reactome:R-HSA-2684901	TAS		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20200515	Reactome
UniProtKB	A1A5B4	ANO9	involved_in	GO:0034220	Reactome:R-HSA-983712	TAS		P	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20210310	Reactome
UniProtKB	A5X5Y0	HTR3E	enables	GO:0022850	PMID:17392525	IDA		F	5-hydroxytryptamine receptor 3E	HTR3E	protein	taxon:9606	20130210	CACAO
UniProtKB	A6NJY1	SLC9B1P1	enables	GO:0015299	GO_REF:0000002	IEA	InterPro:IPR006153	FPutative SLC9B1-like protein SLC9B1P1	SLC9B1P1	protein	taxon:9606	20210612	InterPro
```

### action = gff3

```
$ wget -q -O - "http://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.chr.gff3.gz" | gunzip -c | java -jar dist/goutils.jar --action gff3 -A 'GO:0005216'   | head
WARNING	2021-10-21 16:19:29	AsciiLineReader	Creating an indexable source for an AsciiFeatureCodec using a stream that is neither a PositionalBufferedStream nor a BlockCompressedInputStream
##gff-version 3.1.25
1	ensembl_havana	gene	1215816	1227409	.	+	.	ID=gene%3AENSG00000162572;Name=SCNN1D;biotype=protein_coding;description=sodium channel%2C non-voltage-gated 1%2C delta subunit %5BSource%3AHGNC Symbol%3BAcc%3A10601%5D;gene_id=ENSG00000162572;logic_name=ensembl_havana_gene;version=15
1	ensembl_havana	gene	1950780	1962192	.	+	.	ID=gene%3AENSG00000187730;Name=GABRD;biotype=protein_coding;description=gamma-aminobutyric acid %28GABA%29 A receptor%2C delta %5BSource%3AHGNC Symbol%3BAcc%3A4084%5D;gene_id=ENSG00000187730;logic_name=ensembl_havana_gene;version=6
1	ensembl_havana	gene	6051526	6161253	.	+	.	ID=gene%3AENSG00000069424;Name=KCNAB2;biotype=protein_coding;description=potassium voltage-gated channel%2C shaker-related subfamily%2C beta member 2 %5BSource%3AHGNC Symbol%3BAcc%3A6229%5D;gene_id=ENSG00000069424;logic_name=ensembl_havana_gene;version=10
1	ensembl_havana	gene	11866207	11903201	.	+	.ID=gene%3AENSG00000011021;Name=CLCN6;biotype=protein_coding;description=chloride channel%2C voltage-sensitive 6 %5BSource%3AHGNC Symbol%3BAcc%3A2024%5D;gene_id=ENSG00000011021;logic_name=ensembl_havana_gene;version=17
1	ensembl_havana	gene	13801445	13840543	.	-	.ID=gene%3AENSG00000162494;Name=LRRC38;biotype=protein_coding;description=leucine rich repeat containing 38 %5BSource%3AHGNC Symbol%3BAcc%3A27005%5D;gene_id=ENSG00000162494;logic_name=ensembl_havana_gene;version=5
1	ensembl_havana	gene	16345370	16360545	.	+	.ID=gene%3AENSG00000186510;Name=CLCNKA;biotype=protein_coding;description=chloride channel%2C voltage-sensitive Ka %5BSource%3AHGNC Symbol%3BAcc%3A2026%5D;gene_id=ENSG00000186510;logic_name=ensembl_havana_gene;version=7
1	ensembl_havana	gene	16370272	16383803	.	+	.ID=gene%3AENSG00000184908;Name=CLCNKB;biotype=protein_coding;description=chloride channel%2C voltage-sensitive Kb %5BSource%3AHGNC Symbol%3BAcc%3A2027%5D;gene_id=ENSG00000184908;logic_name=ensembl_havana_gene;version=13
1	ensembl_havana	gene	25071848	25170815	.	+	.ID=gene%3AENSG00000169504;Name=CLIC4;biotype=protein_coding;description=chloride intracellular channel 4 %5BSource%3AHGNC Symbol%3BAcc%3A13518%5D;gene_id=ENSG00000169504;logic_name=ensembl_havana_gene;version=10
1	ensembl_havana	gene	26517052	26529459	.	+	.ID=gene%3AENSG00000188782;Name=CATSPER4;biotype=protein_coding;description=cation channel%2C sperm associated 4 %5BSource%3AHGNC Symbol%3BAcc%3A23220%5D;gene_id=ENSG00000188782;logic_name=ensembl_havana_gene;version=3
```


