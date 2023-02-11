# SamplesRDF

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Digests a  database of samples from a set of recfiles


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samplesrdf  [options] Files

Usage: samplesrdf [options] Files
  Options:
    -doid, --doid
      OWL-formatted Disease ontology. https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.owl
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hpo, --hpo
      OWL-formatted HPO ontology 
      https://github.com/obophenotype/human-phenotype-ontology/blob/master/hp.owl 
  * -o, --out
      Base filename for output files
    -population, --population, --pop, --snomed
      OWL-formatted SNOWMED-ETHNIC-GROUP 
      https://bioportal.bioontology.org/ontologies/SNOMED-Ethnic-Group 
    --version
      print version and exit

```


## Keywords

 * samples
 * database
 * rdf
 * ontology



## Creation Date

20230201

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samplesrdf/SamplesRDF.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samplesrdf/SamplesRDF.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samplesrdf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

A tool formatting samples definition.

# Input Format

Input is a Recfile ( https://en.wikipedia.org/wiki/Recfiles )  in the form:

```
%rec: type1

key : value
key : value
key : value
key : value

%rec: type2

key : value
key : value
key : value
key : value

key : value
(...)
```

the key/value delimiter is `:` 
multiple records are separated by one or more spaces.
lines starting with "#" are ignored.
keys are case insensible.

A rec file must contain a class declaration wich introduce the classes for the following object:
```
%rec: Sample
```
or
```
%rec: Group
```

The type is either Sample or Group. 
A Sample defines one Sample.
A group is a set of samples.



**id:** : each record requires the key `id`. The type is either Sample or Group. It defines a unique identifier in the database. Id are converted to uppercase

**label:** : a title for this record

**comment:** or '**description:** : a description for the record

** hpo: ** : or '**hp:** is a identifier in the human phenotype ontology database. If the `hpo:` is defined in a `type:Group`, the phenotype will be propaged to all the samples in that group.
The following syntaxes are allowed:
```
hp: HP:0011712 everything after the first word is ignored
hp: 0011712 words
hpo: 0011712 words
```


** doid:** is a identifier in the human disease ontology database. If the `doid:` is defined in a `type:Group`, the phenotype will be propaged to all the samples in that group.
The following syntaxes are allowed:
```
doid: DOID:0011712 everything after the first word is ignored
doid: 0011712 words
doid: 0011712 words
```



## Samples only:

**family:** a family name for this a `Sample`

**father:** the id of the father 

**sex:** the sex of the sample 'male' or 'female'

**birth:** the  year of birth

**alias:** another name for this sample


## group:

**group: id ** or **extends:id** : the current group extends another group defined by 'id'


**pop** or **population** is an identifier in the SNOMED population ontology. https://bioportal.bioontology.org/ontologies/SNOMED-Ethnic-Grou . Value is either the label or the SNOMEDID in the OWL ontology
```
type: sample
id: John
pop: Danes
pop: S-61240
```


 
# Output Format

 * (prefix)_errors.txt : errors found
 * (prefix)_model.rdf : output formatted as RDF/XML
 * (prefix)_doid.txt: output each item in Disease Ontology. Each line is a term in DO listing each sample with this disease.
 * (prefix)_pedigree.txt : output samples formatted as a pedigree file. The links / sex are not checked
 * (prefix)_samples.txt : output samples with their diseases, phenotypes...


# Example

```
$ cat input.recfile 

type:sample
id: azd
father: x1
sex: male
doid: DOID:0050451 Hello
hpo: HP:0011712 World
pop: S-61020

type:group
id:g1
sample: azd
sample: x1


 java -jar dist/jvarkit.jar samplesrdf --population  pop.owl --doid doid.owl --hpo hp.owl -o JETER input.recfile
 
 $ more JETER_doid.txt 
http://purl.obolibrary.org/obo/DOID_0050451	DOID:0050451	Brugada syndrome	AZD
http://purl.obolibrary.org/obo/DOID_0110221	DOID:0110221	Brugada syndrome 4	AZD
http://purl.obolibrary.org/obo/DOID_0110222	DOID:0110222	Brugada syndrome 5	AZD
http://purl.obolibrary.org/obo/DOID_0110220	DOID:0110220	Brugada syndrome 3	AZD
http://purl.obolibrary.org/obo/DOID_0110225	DOID:0110225	Brugada syndrome 8	AZD
http://purl.obolibrary.org/obo/DOID_0110226	DOID:0110226	Brugada syndrome 9	AZD
http://purl.obolibrary.org/obo/DOID_0110223	DOID:0110223	Brugada syndrome 6	AZD
http://purl.obolibrary.org/obo/DOID_0110224	DOID:0110224	Brugada syndrome 7	AZD
http://purl.obolibrary.org/obo/DOID_0110219	DOID:0110219	Brugada syndrome 2	AZD
http://purl.obolibrary.org/obo/DOID_0110218	DOID:0110218	Brugada syndrome 1	AZD


$ more JETER_hpo.txt 
http://purl.obolibrary.org/obo/HP_0011712	HP:0011712	Right bundle branch block	AZD

$ more JETER_groups.txt 
G1	G1		2	1	0	AZD;X1


$ more JETER_samples.txt
X1	X1		.	.	.	.	.	.	.	.	G1
AZD	AZD		.	.	X1	.	male	Right bundle branch block	Brugada syndrome 8; Brugada syndrome 7; Brugada syndrome 6; Brugada syndrome 5; Brugada syndrome 1; B
rugada syndrome; Brugada syndrome 9; Brugada syndrome 2; Brugada syndrome 4; Brugada syndrome 3	Armenians	G1

$ more JETER_pedigree.txt 
X1	X1	0	0	0
AZD	AZD	X1	0	1

$ more JETER_errors.txt 
<https://umr1087.univ-nantes.fr/db/X1> rdf:type <https://umr1087.univ-nantes.fr/Sample> was used but not declared
X1 declared as father of AZD. But sex is not male

$ more JETER_model.rdf 
<rdf:RDF
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:u="https://umr1087.univ-nantes.fr/" > 
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/G1">
    <u:contains-sample rdf:resource="https://umr1087.univ-nantes.fr/db/AZD"/>
    <u:id>G1</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Group"/>
    <u:contains-sample rdf:resource="https://umr1087.univ-nantes.fr/db/X1"/>
  </rdf:Description>
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/X1">
    <u:id>X1</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Sample"/>
  </rdf:Description>
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/AZD">
    <u:hpo rdf:resource="http://purl.obolibrary.org/obo/HP_0011712"/>
    <u:doid rdf:resource="http://purl.obolibrary.org/obo/DOID_0050451"/>
    <u:sex rdf:resource="https://umr1087.univ-nantes.fr/Male"/>
    <u:father rdf:resource="https://umr1087.univ-nantes.fr/db/X1"/>
    <u:id>AZD</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Sample"/>
    <u:population rdf:resource="http://purl.bioontology.org/ontology/SNOMED-Ethnic-Group#113169009"/>
  </rdf:Description>
</rdf:RDF>
```


