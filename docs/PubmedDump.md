# PubmedDump

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Dump XML results from pubmed/Eutils


## Usage

```
Usage: pubmeddump [options] Files
  Options:
    -e, --email
      optional user email
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ncbi-api-key
      NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 
      .If undefined, it will try to get in that order:  1) environment 
      variable ${NCBI_API_KEY} ;  2) the jvm property "ncbi.api.key" ;	3) A 
      java property file ${HOME}/.ncbi.properties and key api_key
    -o, --output
      Output file. Optional . Default: stdout
    -skip, --skip
      [20180302]  Optional set of elements names to be ignored in the output. 
      Spaces or comma separated. .eg: 'AuthorList PubmedData '
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * ncbi
 * pubmed
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/270498](https://www.biostars.org/p/270498)
 * [https://www.biostars.org/p/365479](https://www.biostars.org/p/365479)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmeddump
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedDump.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedDump.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedDumpTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedDumpTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmeddump** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$  java -jar dist/pubmeddump.jar "Lindenbaum P" | grep ArticleTitle
			
	<ArticleTitle>NGS library preparation may generate artifactual integration sites of AAV vectors.</ArticleTitle>
    <ArticleTitle>Mutations in FAM111B cause hereditary fibrosing poikiloderma with tendon contracture, myopathy, and pulmonary fibrosis.</ArticleTitle>
    <ArticleTitle>[The Spanish Association of Surgeon's audited teaching programme for rectal cancer. Results after six years].</ArticleTitle>
    <ArticleTitle>Common variants at SCN5A-SCN10A and HEY2 are associated with Brugada syndrome, a rare disease with high risk of sudden cardiac death.</ArticleTitle>
    <ArticleTitle>The 3rd DBCLS BioHackathon: improving life science data integration with Semantic Web technologies.</ArticleTitle>
    <ArticleTitle>Mass spectrometry-based identification of native cardiac Nav1.5 channel alpha subunit phosphorylation sites.</ArticleTitle>
    <ArticleTitle>BioStar: an online question &amp; answer resource for the bioinformatics community.</ArticleTitle>
    <ArticleTitle>Knime4Bio: a set of custom nodes for the interpretation of next-generation sequencing data with KNIME.</ArticleTitle>
    <ArticleTitle>Truncating mutations in the last exon of NOTCH2 cause a rare skeletal disorder with osteoporosis.</ArticleTitle>
    <ArticleTitle>The Gene Wiki: community intelligence applied to human gene annotation.</ArticleTitle>
    <ArticleTitle>Robust physical methods that enrich genomic regions identical by descent for linkage studies: confirmation of a locus for osteogenesis imperfecta.</ArticleTitle>
    <ArticleTitle>Association of autism with polymorphisms in the paired-like homeodomain transcription factor 1 (PITX1) on chromosome 5q31: a candidate gene analysis.</ArticleTitle>
    <ArticleTitle>Haplotypes in the gene encoding protein kinase c-beta (PRKCB1) on chromosome 16 are associated with autism.</ArticleTitle>
    <ArticleTitle>RoXaN, a novel cellular protein containing TPR, LD, and zinc finger motifs, forms a ternary complex with eukaryotic initiation factor 4G and rotavirus NSP3.</ArticleTitle>
    <ArticleTitle>CloneIt: finding cloning strategies, in-frame deletions and frameshifts.</ArticleTitle>
    <ArticleTitle>In vivo and in vitro phosphorylation of rotavirus NSP5 correlates with its localization in viroplasms.</ArticleTitle>
```

## See also

 * https://gist.github.com/lindenb/6bfb49fd8bc3dd27d99f
 * https://gist.github.com/lindenb/5d7773a93d8c2b0edbd4c01bf8834919

