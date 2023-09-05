# RDFCombine

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Substract/Add RDF models


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar rdfcombine  [options] Files

Usage: rdfcombine [options] Files
  Options:
    --base
      (URI) xml:base when reading rdf model from stdin or writing model
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --keep-rdfType
      ignore rdf:type in minus model. (do not remove rdf:type in source)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --output-format
      Write a serialized represention of a model in a specified language. 
      Predefined values are RDF/XML, RDF/XML-ABBREV, N-TRIPLE, TURTLE, (and 
      TTL) and N3
      Default: RDF/XML-ABBREV
    --version
      print version and exit

```


## Creation Date

20230903

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/rdfcombine/RDFCombine.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/rdfcombine/RDFCombine.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **rdfcombine** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a list of URI, or a file with the '.list' suffix containing the URI

URI are processed sequentially.

URI can be prefixed with '+' or '-' . A minus sign means that the model will be substracted to the current one.

WARNING: blank nodes are always not equal

## Example

```
$ java -jar dist/jvarkit.jar rdfcombine ${HOME}/file.rdf   -${HOME}/file.rdf 
[INFO][RDFCombine]number of statements after adding file.rdf) = 1663
[INFO][RDFCombine]number of statements after substracting -file.rdf) = 0
<rdf:RDF
    xmlns:rel="http://purl.org/vocab/relationship/"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:bio="http://purl.org/vocab/bio/0.1/"
    xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
    xmlns:foaf="http://xmlns.com/foaf/0.1/">
</rdf:RDF>


$ java -jar dist/jvarkit.jar rdfcombine ${HOME}/file2.rdf  ${HOME}/file2.rdf 


java -jar dist/jvarkit.jar  rdfcombine "https://raw.githubusercontent.com/BruceMWhealton/Gedcom-RDF/master/Disney.rdf" "https://raw.githubusercontent.com/BruceMWhealton/Gedcom-RDF/master/Thomas.rdf"

```


