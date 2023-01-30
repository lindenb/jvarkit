# OboUtils

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

OBO Ontology Utils.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar oboutils  [options] Files

Usage: oboutils [options] Files
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
      Possible Values: [dump_table]
    --exclude-accession, --exclude
      User Go Terms to be EXCLUDED accession numbers or name.eg
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --inverse
      inverse the result
      Default: false
    -obo, --obo, --ontology
      Ontology in OBO format
      Default: https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * obo
 * ontology



## Creation Date

20230105

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/obo/OboUtils.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/obo/OboUtils.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **oboutils** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/oboutils.jar --obo "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo" -A "HP:0001626" | head
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




