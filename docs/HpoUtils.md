# HpoUtils

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Human Phenotype Ontology Utils.


## Usage

```
Usage: java -jar dist/hpoutils.jar  [options] Files
Usage: hpoutils [options] Files
  Options:
    -A, --accession
      User Go Terms accession numbers or name.eg HP:0001626 Abnormality of the 
      cardiovascular system
      Default: []
    -af, --accession-file
      File containing accession numbers. One per line.
    -action, --action
      What shoud I do ? default is dump as table.
      Default: dump_table
      Possible Values: [dump_table, genes, gff3]
    --exclude-accession, --exclude
      User Go Terms to be EXCLUDED accession numbers or name.eg
      Default: []
    -g2p, --g2p
      HPO gene to phenotype file. 
      eg://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt ( Format: 
      entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-ID<tab> )
    -g, --gff, --gff3
      GFF3 file for action=genes or action=gff3.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hpo, --hpo, --obo
      HPO ontology in OBO format
      Default: https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
    -i, --inverse
      inverse the result
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * phenotype
 * ontology
 * hpo
 * hpoa


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew hpoutils
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20230105

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/hpo/HpoUtils.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/hpo/HpoUtils.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **hpoutils** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/hpoutils.jar -A "HP:0001626" | head
#ACN	NAME	DEFINITION
HP:0001195	Single umbilical artery	Single umbilical artery (SUA) is the absence of one of the two umbilical arteries surrounding the fetal bladder and in the fetal umbilical cord.
HP:0001136	Retinal arteriolar tortuosity	The presence of an increased number of twists and turns of the retinal arterioles.
HP:0025188	Retinal vasculitis	Inflammation of retinal blood vessels as manifested by perivascular sheathing or cuffing, vascular leakage and/or occlusion.
HP:0025169	Left ventricular systolic dysfunction	Abnormality of left ventricular contraction, often defined operationally as an ejection fraction of less than 40 percent.
HP:0025168	Left ventricular diastolic dysfunction	Abnormal function of the left ventricule during left ventricular relaxation and filling.
HP:0410173	Increased circulating troponin I concentration	An increased concentration of tropnin I in the blood, which is a cardiac regulatory protein that controls the calcium mediated interaction between actin and myosin. Raised cardiac troponin concentrations are now accepted as the standard biochemical marker for the diagnosis of myocardial infarction.
HP:0410174	Increased circulating troponin T concentration	An increased concentration of tropnin T in the blood, which is a cardiac regulatory protein that controls the calcium mediated interaction between actin and myosin. Raised cardiac troponin concentrations are now accepted as the standard biochemical marker for the diagnosis of myocardial infarction.
HP:0410267	Intestinal hemangioma	A hemangioma, a benign tumor of the vascular endothelial cells, located in the intestines, which includes the bowel.
HP:0410268	Spleen hemangioma	A hemangioma, a benign tumor of the vascular endothelial cells, that is located in the spleen.
```

