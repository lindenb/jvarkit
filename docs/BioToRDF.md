# BioToRDF

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Build a RDF database for human from misc sources


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bio2rdf  [options] Files

Usage: bio2rdf [options] Files
  Options:
    --genes
      Limit to those genes names , separated with comma (for debugging)
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -D
      parameters. -Dkey1=value1  -Dkey2=value2 ...
      Syntax: -Dkey=value
      Default: {GO_OWL=http://purl.obolibrary.org/obo/go.owl, HUMAN_HPO_OWL=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/hp.owl, HPO_PHENOTYPE_TO_GENE=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/phenotype_to_genes.txt, NCBI_GENE_INFO=https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz, GFF3_GRCH38=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.annotation.gff3.gz, GFF3_GRCH37=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/GRCh37_mapping/gencode.v{GENCODE_RELEASE}lift37.annotation.gff3.gz, GENCODE_RELEASE=43, MONDO_OWL=http://purl.obolibrary.org/obo/mondo.owl, NCBI_GENE_GO=https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz}

```


## Keywords

 * rdf
 * ontology
 * sparql



## Creation Date

20220427

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bio2rdf/BioToRDF.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bio2rdf/BioToRDF.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bio2rdf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# motivation

build  RDF file to build a RDF database.
 
 
# usage
 
```
$ java -jar dist/jvarkit.jar bio2rdf | gzip > bio2rdf.rdf.gz
```

## Example queries:

### Example

find all the descendant of `GO:0045823`

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>


SELECT ?subClass ?label ?desc WHERE { 
        ?subClass rdfs:subClassOf* <http://purl.obolibrary.org/obo/GO_0045823> . 
        ?subClass dc:title ?label . 
        ?subClass dc:description ?desc . 
    }
```

```
$ arq --data=/home/me/bio2rdf.rdf.gz --query query.sparql
```

| subClass                                    | label        | desc                                                                                                                                    |
|--- |-- |-- |
| <http://purl.obolibrary.org/obo/GO_0045823> | "GO:0045823" | "positive regulation of heart contraction"                                                                                              |
| <http://purl.obolibrary.org/obo/GO_0001989> | "GO:0001989" | "positive regulation of the force of heart contraction involved in baroreceptor response to decreased systemic arterial blood pressure" |
| <http://purl.obolibrary.org/obo/GO_0010460> | "GO:0010460" | "positive regulation of heart rate"                                                                                                     |
| <http://purl.obolibrary.org/obo/GO_0001996> | "GO:0001996" | "positive regulation of heart rate by epinephrine-norepinephrine"                                                                       |
| <http://purl.obolibrary.org/obo/GO_0086024> | "GO:0086024" | "adenylate cyclase-activating adrenergic receptor signaling pathway involved in positive regulation of heart rate"                      |
| <http://purl.obolibrary.org/obo/GO_0003065> | "GO:0003065" | "positive regulation of heart rate by epinephrine"                                                                                      |
| <http://purl.obolibrary.org/obo/GO_0003112> | "GO:0003112" | "positive regulation of heart rate by neuronal epinephrine"                                                                             |
| <http://purl.obolibrary.org/obo/GO_0003111> | "GO:0003111" | "positive regulation of heart rate by circulating epinephrine"                                                                          |
| <http://purl.obolibrary.org/obo/GO_0003066> | "GO:0003066" | "positive regulation of heart rate by norepinephrine"                                                                                   |
| <http://purl.obolibrary.org/obo/GO_0003114> | "GO:0003114" | "positive regulation of heart rate by circulating norepinephrine"                                                                       |
| <http://purl.obolibrary.org/obo/GO_0003113> | "GO:0003113" | "positive regulation of heart rate by neuronal norepinephrine"                                                                          |
| <http://purl.obolibrary.org/obo/GO_0001988> | "GO:0001988" | "positive regulation of heart rate involved in baroreceptor response to decreased systemic arterial blood pressure"                     |
| <http://purl.obolibrary.org/obo/GO_0060452> | "GO:0060452" | "positive regulation of cardiac muscle contraction"                                                                                     |
| <http://purl.obolibrary.org/obo/GO_0003099> | "GO:0003099" | "positive regulation of the force of heart contraction by chemical signal"                                                              |
| <http://purl.obolibrary.org/obo/GO_0003061> | "GO:0003061" | "positive regulation of the force of heart contraction by norepinephrine"                                                               |
| <http://purl.obolibrary.org/obo/GO_0003110> | "GO:0003110" | "positive regulation of the force of heart contraction by neuronal norepinephrine"                                                      |
| <http://purl.obolibrary.org/obo/GO_0003109> | "GO:0003109" | "positive regulation of the force of heart contraction by circulating norepinephrine"                                                   |
| <http://purl.obolibrary.org/obo/GO_0003059> | "GO:0003059" | "positive regulation of the force of heart contraction by epinephrine"                                                                  |
| <http://purl.obolibrary.org/obo/GO_0003087> | "GO:0003087" | "positive regulation of the force of heart contraction by neuronal epinephrine"                                                         |
| <http://purl.obolibrary.org/obo/GO_0003088> | "GO:0003088" | "positive regulation of the force of heart contraction by circulating epinephrine"                                                      |
| <http://purl.obolibrary.org/obo/GO_0001997> | "GO:0001997" | "positive regulation of the force of heart contraction by epinephrine-norepinephrine"                                                   |
| <http://purl.obolibrary.org/obo/GO_0003090> | "GO:0003090" | "positive regulation of the force of heart contraction by neuronal epinephrine-norepinephrine"                                          |
| <http://purl.obolibrary.org/obo/GO_0003089> | "GO:0003089" | "positive regulation of the force of heart contraction by circulating epinephrine-norepinephrine"                                       |

### find all the phenotypes containing the word 'arrhythmia'

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX bio: <https://umr1087.univ-nantes.fr/bio2rdf/>

SELECT  ?label ?desc WHERE { 
        ?disease a bio:Phenotype .
        ?disease dc:title ?label . 
        ?disease dc:description ?desc . 
	FILTER (CONTAINS(LCASE(?desc),"arrhythmia")) .
    }
```

```
------------------------------------------------
| label        | desc                          |
================================================
| "HP:0002521" | "Hypsarrhythmia"              |
| "HP:0011215" | "Hemihypsarrhythmia"          |
| "HP:0004308" | "Ventricular arrhythmia"      |
| "HP:0001692" | "Atrial arrhythmia"           |
| "HP:0005115" | "Supraventricular arrhythmia" |
| "HP:0011675" | "Arrhythmia"                  |
------------------------------------------------
```


