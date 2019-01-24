# PubmedOrcidGraph

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Creates a graph from Pubmed and Authors' Orcid identifiers


## Usage

```
Usage: pubmedorcidgraph [options] Files
  Options:
    -links, --alllinks
      By default, we display only one link between two authors. Using this 
      option will show all the links (publications)
      Default: false
  * -D, --berkeydb
      BerkeleyDB tmpDir
    -ea, --edgeattributes
      Do not show edge attributes (smaller files with less informations)
      Default: false
    -E, --errors
      Dump strange orcids (e.g: same orcird but different forename) . Default: 
      stderr 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -d, --maxdepth
      Max graph depth
      Default: 2
    --ncbi-api-key
      NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 
      .If undefined, it will try to get in that order:  1) environment 
      variable ${NCBI_API_KEY} ;  2) the jvm property "ncbi.api.key" ;	3) A 
      java property file ${HOME}/.ncbi.properties and key api_key
    -orcid, --orcid
      Input is a set of orcids identifiers
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * pubmed
 * ncbi
 * orcid


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmedorcidgraph
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedOrcidGraph.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedOrcidGraph.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedorcidgraph** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## About Orcid
 ORCID  (http://orcid.org) provides a persistent digital identifier that distinguishes you from every other researcher and, through integration in key research workflows such as manuscript and grant submission, supports automated linkages between you and your professional activities ensuring that your work is recognized.

You can download the papers containing some orcid Identifiers using the entrez query

http://www.ncbi.nlm.nih.gov/pubmed/?term=orcid[AUID]

I've used one of my tools pubmeddump to download the articles as XML and I wrote PubmedOrcidGraph to extract the author's orcid.
The output is a GEXF file for gephi.

## Example

using pubmed efetch output

```
java -jar dist/pubmeddump.jar --skip "MeshHeadingList ChemicalList GrantList InvestigatorList CommentsCorrectionsList ISSN DateRevised AffiliationInfo Language PublicationTypeList  ArticleDate PubmedData Abstract MedlineJournalInfo CoiStatement KeywordList Pagination ELocationID "   "orcid[AUID]" |\
java -jar dist/pubmedorcidgraph.jar -D BDB 
```

using orcid identifiers:

java -jar dist/pubmedorcidgraph.jar -D BDB --orcid 0000-0001-7751-2280 0000-0003-0677-5627 0000-0003-3628-1548 0000-0003-4530-6655 0000-0001-8007-5931 


<img src="https://pbs.twimg.com/media/Ci-h0MJWUAAvJjw.jpg"/>

## See also</h:h3>

 * http://plindenbaum.blogspot.fr/2016/05/playing-with-orcidorg-ncbipubmed-graph.html

