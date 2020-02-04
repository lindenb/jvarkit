# Pubmed404

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Test if URL in the pubmed abstracts are reacheable.


## Usage

```
Usage: pubmed404 [options] Files
  Options:
    -c, --collapse
      Only one URL per article. Print the '200/OK' first.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -t, --timeout
      timeout in seconds
      Default: 5
    --version
      print version and exit

```


## Keywords

 * pubmed
 * url


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmed404
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20181210

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/Pubmed404.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/Pubmed404.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmed404** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/pubmeddump.jar 'bioinformatics 2001' 2> /dev/null |\
	java -jar dist/pubmed404.jar  2> /dev/null 

#PMID	TITLE	YEAR	URL	Status
29520589	Expression of Colocasia esculenta tuber agglutinin in Indian mustard provides resistance against Lipaphis erysimi and the expressed protein is non-allergenic.2018	http://www.fao.org/docrep/007/y0820e/y0820e00.HTM	200
29520589	Expression of Colocasia esculenta tuber agglutinin in Indian mustard provides resistance against Lipaphis erysimi and the expressed protein is non-allergenic.2018	http://www.icmr.nic.in/guide/Guidelines%20for%20Genetically%20Engineered%20Plants.pdf	-1
28482857	Horizontal gene transfer is not a hallmark of the human genome.	2017	https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0607-3	200
27899642	The UCSC Genome Browser database: 2017 update.	2017	http://genome.ucsc.edu/	200
27797935	High hospital research participation and improved colorectal cancer survival outcomes: a population-based study.	2017	http://www.bmj.com/company/products-services/rights-and-licensing/	403
25505092	NMRFAM-SPARKY: enhanced software for biomolecular NMR spectroscopy.	2015	http://pine.nmrfam.wisc.edu/download_packages.html	200
25505092	NMRFAM-SPARKY: enhanced software for biomolecular NMR spectroscopy.	2015	http://www.nmrfam.wisc.edu/nmrfam-sparky-distribution.htm	200
25428374	The UCSC Genome Browser database: 2015 update.	2015	http://genome.ucsc.edu	200
26356339	A Simple but Powerful Heuristic Method for Accelerating k-Means Clustering of Large-Scale Data in Life Science.	null	http://mlab.cb.k.u-tokyo.ac.jp/~ichikawa/boostKCP/	200
24794704	Usefulness of the Shock Index as a secondary triage tool.	2015	http://group.bmj.com/group/rights-licensing/permissions	403
24225322	Progenetix: 12 years of oncogenomic data curation.	2014	http://www.progenetix.org	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://globin.bx.psu.edu/hbvar	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://www.findbase.org	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://www.lovd.nl	200
23564938	DAMBE5: a comprehensive software package for data analysis in molecular biology and evolution.	2013	http://dambe.bio.uottawa.ca	200
22689647	SIFT web server: predicting effects of amino acid substitutions on proteins.	2012	http://sift-dna.org	200
22600740	Cyber-T web server: differential analysis of high-throughput data.	2012	http://cybert.ics.uci.edu/	200
21742331	An open source lower limb model: Hip joint validation.	2011	https://simtk.org/home/low_limb_london	200
21593132	Java bioinformatics analysis web services for multiple sequence alignment--JABAWS:MSA.	2011	http://www.compbio.dundee.ac.uk/jabaws	200
20228129	DensiTree: making sense of sets of phylogenetic trees.	2010	http://compevol.auckland.ac.nz/software/DensiTree/	404
19380317	CELLULAR OPEN RESOURCE (COR): current status and future directions.	2009	http://www.cellml.org/specifications/	200
18948284	OperonDB: a comprehensive database of predicted operons in microbial genomes.	2009	http://operondb.cbcb.umd.edu	200
18368364	Simulator for neural networks and action potentials.	2007	http://snnap.uth.tmc.edu	-1
18367465	An improved general amino acid replacement matrix.	2008	http://atgc.lirmm.fr/LG	404
18238804	Interoperability with Moby 1.0--it's better than sharing your toothbrush!	2008	http://www.biomoby.org/	200
18174178	PRALINETM: a strategy for improved multiple alignment of transmembrane proteins.	2008	http://www.ibi.vu.nl/programs/pralinewww	200
17221864	HbVar database of human hemoglobin variants and thalassemia mutations: 2007 update.	2007	http://globin.bx.psu.edu/hbvar	200
17221864	HbVar database of human hemoglobin variants and thalassemia mutations: 2007 update.	2007	http://www.goldenhelix.org/xprbase	403
(...)
```
