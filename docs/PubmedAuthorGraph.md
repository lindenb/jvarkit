# PubmedAuthorGraph

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates a graph from Pubmed and Authors


## Usage

```
Usage: pubmedauthorgraph [options] Files
  Options:
    -rc, --article-color
      viz:Color for the Articles. A named color ('red', 'blue'...) use the 
      syntax 'rgb(int,int,int)'.
    -uc, --author-color
      viz:Color for the Authors. A named color ('red', 'blue'...) use the 
      syntax 'rgb(int,int,int)'.
  * -D, --berkeydb
      BerkeleyDB tmpDir
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --initals
      use author's initials to build the author-identifier. In the old pubmed 
      record, the forename is not available.
      Default: false
    --ncbi-api-key
      NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 
      .If undefined, it will try to get in that order:  1) environment 
      variable ${NCBI_API_KEY} ;  2) the jvm property "ncbi.api.key" ;	3) A 
      java property file ${HOME}/.ncbi.properties and key api_key
    -o, --output
      Output file. Optional . Default: stdout
    -r, --scale-articles
      Scale articles in function of the number of authors
      Default: false
    -u, --scale-authors
      Scale authors in function of the number of articles
      Default: false
    -S, -singleton, --singleton
      Remove singleton authors (authors linked to only one paper)
      Default: false
    --version
      print version and exit

```


## Keywords

 * pubmed
 * ncbi
 * graph


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmedauthorgraph
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedAuthorGraph.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedAuthorGraph.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedauthorgraph** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

builds a graph with XML pubmed. Nodes are 'Article' and  'Authors'. Edges are authorships.

## Example

```
java -jar dist/pubmeddump.jar 'SCN5A redon' |\
java -jar dist/pubmedauthorgraph.jar -D BDB -singleton   > out.gexf
```

## Screenshots

https://twitter.com/yokofakun/status/1034107797439504384

![https://twitter.com/yokofakun/status/1034107797439504384](https://pbs.twimg.com/media/DlnjTj1W4AIBB6B.jpg)

https://twitter.com/yokofakun/status/1034397660189523968

![https://twitter.com/yokofakun/status/1034397660189523968](https://pbs.twimg.com/media/DlrqXqvX4AE6r27.jpg)



