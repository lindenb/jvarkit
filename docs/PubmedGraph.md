# PubmedGraph

Creates a Gephi-gexf graph of references-cotes for a given PMID


## Usage

```
Usage: pubmedgraph [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --ncbi-api-key
      NCBI API Key. If undefined. Will try to read in that order: 1) A java 
      XML property file ${HOME}/.ncbi.properties 2) the jvm property 
      "ncbi.api.key" 3) environment variable NCBI_API_KEY
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -b
      disable backward (referenced-in)
      Default: false
    -d
      max-depth
      Default: 3
    -f
      disable forward (cited-in)
      Default: false

```


## Keywords

 * pubmed
 * xml
 * graph


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) and avoid OpenJdk, use the java from Oracle. Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make pubmedgraph
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGraph.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGraph.java)


<details>
<summary>Git History</summary>

```
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Tue Apr 4 17:09:36 2017 +0200 ; vcfgnomad ; https://github.com/lindenb/jvarkit/commit/eac33a01731eaffbdc401ec5fd917fe345b4a181
Thu Jul 28 09:48:29 2016 +0200 ; NCBI moved API to https ; https://github.com/lindenb/jvarkit/commit/d207e023a06d2ae7afd2e05d2f1369b8a713974b
Tue Jun 16 17:40:00 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/f394e306c86bd45240d165c69748acf44f0b38ec
Mon Jun 15 21:35:36 2015 +0200 ; pubmed graph ; https://github.com/lindenb/jvarkit/commit/3fd6216961e9aca27a1013a7b2aab63afe819d52
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedgraph** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


