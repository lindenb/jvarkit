# EVADumpFiles

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Dump files locations from European Variation Archive


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar evadumpfiles  [options] Files

Usage: evadumpfiles [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --taxons, --taxon
      define one or more taxon filters:syntax 'taxonCode assemblyCode' or 
      'taxonCode *' or 'taxonCode' .eg: 'cfamiliaris' , 'cfamiliaris * ', 
      'cfamiliaris 31'. Multiple are comma/semicolon separated. See output of 
      https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/species/list/ .
      Default: []
    --version
      print version and exit

```


## Keywords

 * eva
 * ebi
 * snp
 * variant



## Creation Date

20230314

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/eva/EVADumpFiles.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/eva/EVADumpFiles.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **evadumpfiles** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example
```
$ java -jar dist/jvarkit.jar evadumpfiles --taxons 'cfamiliaris 31' 2> /dev/null |  verticalize  

>>> 2
$1	studyId	PRJEB56315
$2	studyName	Idiopathic Epilepsy in the Dutch partridge dog
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB56315/variants_pass.vcf.gz
<<< 2

>>> 3
$1	studyId	PRJEB56211
$2	studyName	Performance of Variant Pathogenicity Prediction Methods
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB56211/VPP_Project_Variants.vcf.gz
<<< 3

>>> 4
$1	studyId	PRJEB51024
$2	studyName	Longevity of Cane Corso Italiano purebred dogs
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB51024/cane_corso_f.vcf.gz
```


