# VcfSparql

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Query RDF with Sparql. Very slow. Just a proof of concept


## Usage

```
Usage: vcfsparql [options] Files
  Options:
    --code
      show code
      Default: false
    --format
      output format. One of https://jena.apache.org/documentation/javadoc/arq/org/apache/jena/sparql/resultset/ResultsFormat.html#lookup-java.lang.String-
      Default: text
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -a, --prefixes
      Prepend common RDF prefixes to the expression
      Default: false
    -r, --region
      limit query to this genomic interval. An interval as the following 
      syntax : "chrom:start-end" or "chrom:middle+extend"  or 
      "chrom:start-end+extend" or "chrom:start-end+extend-percent%".A program 
      might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)
    -e, --sparql
      Query as sparql string
    -f, --sparql-file
      Query as sparql file
    --version
      print version and exit

```


## Keywords

 * vcf
 * sparql
 * rdf
 * arq
 * semanticweb


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfsparql
```

The java jar file will be installed in the `dist` directory.


## Creation Date

2019-03-06

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsparql/VcfSparql.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsparql/VcfSparql.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfsparql/VcfSparqlTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfsparql/VcfSparqlTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsparql** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example
 
``` 
$ cat query.sparql


PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>


SELECT distinct *
WHERE {
    ?s ?p ?o .
    }
```

```
$ java -jar dist/vcfsparql.jar -f query.sparql src/test/resources/rotavirus_rf.vcf.gz

---------------------------------------------------------------------------------------------------------------
| s      | p                                                 | o                                              |
===============================================================================================================
| _:b0   | <jvarkit:allele>                                  | "C"                                            |
| _:b0   | <jvarkit:sample>                                  | "S5"                                           |
| _:b1   | <jvarkit:type>                                    | "HOM_VAR"                                      |
| _:b1   | <jvarkit:genotype>                                | _:b0                                           |
| _:b0   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b2   | <jvarkit:allele>                                  | "A"                                            |
| _:b2   | <jvarkit:sample>                                  | "S4"                                           |
| _:b1   | <jvarkit:type>                                    | "HOM_REF"                                      |
| _:b1   | <jvarkit:genotype>                                | _:b2                                           |
| _:b2   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b3   | <jvarkit:allele>                                  | "A"                                            |
| _:b3   | <jvarkit:sample>                                  | "S3"                                           |
| _:b1   | <jvarkit:genotype>                                | _:b3                                           |
| _:b3   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b4   | <jvarkit:allele>                                  | "A"                                            |
| _:b4   | <jvarkit:sample>                                  | "S2"                                           |
| _:b1   | <jvarkit:genotype>                                | _:b4                                           |
| _:b4   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b5   | <jvarkit:allele>                                  | "A"                                            |
| _:b5   | <jvarkit:sample>                                  | "S1"                                           |
(...)
---------------------------------------------------------------------------------------------------------------
```

 
