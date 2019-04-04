# PubmedMap

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Use Pubmed Author's Affiliation to map the authors in the world.


## Usage

```
Usage: pubmedmap [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * pubmed
 * xml
 * gis
 * map


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmedmap
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedMap.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedMap.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedMapTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedMapTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedmap** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation
use the Affiliation field in XML pubmed and try to insert some XML attributes describing the location.

## Example

```
$ curl -s "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=26941271&retmode=xml" |java -jar dist/pubmedmap.jar   |\
grep domain -B 5

                 <Author>
                    <LastName>Lai</LastName>
                    <ForeName>Chih-Cheng</ForeName>
                    <Initials>CC</Initials>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Intensive Care Medicine, Chi Mei Medical Center, Liouying, Tainan, Taiwan.</Affiliation>
--
                    <LastName>Lee</LastName>
                    <ForeName>Meng-Tse Gabriel</ForeName>
                    <Initials>MG</Initials>
                    <Identifier Source="ORCID">http://orcid.org/0000-0002-2648-1522</Identifier>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Emergency Medicine, National Taiwan University Hospital, Taipei, Taiwan.</Affiliation>
--
                    <LastName>Lee</LastName>
                    <ForeName>Shih-Hao</ForeName>
                    <Initials>SH</Initials>
                    <Identifier Source="ORCID">http://orcid.org/0000-0002-2648-1522</Identifier>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Emergency Medicine, National Taiwan University Hospital, Taipei, Taiwan.</Affiliation>
--
                <Author>
                    <LastName>Hsu</LastName>
                    <ForeName>Wan-Ting</ForeName>
                    <Initials>WT</Initials>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Emergency Medicine, National Taiwan University Hospital, Taipei, Taiwan.</Affiliation>
--
                <Author>
                    <LastName>Chang</LastName>
                    <ForeName>Shy-Shin</ForeName>
                    <Initials>SS</Initials>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Family Medicine, Chang Gung Memorial Hospital, Linkou, Taiwan Graduate Institute of Clinical Medical Sciences, College of Medicine, Chang Gung University, Taoyuan, Taiwan.</Affiliation>
--
                <Author>
                    <LastName>Chen</LastName>
                    <ForeName>Shyr-Chyr</ForeName>
                    <Initials>SC</Initials>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Emergency Medicine, National Taiwan University Hospital, Taipei, Taiwan.</Affiliation>
--
                    <LastName>Lee</LastName>
                    <ForeName>Chien-Chang</ForeName>
                    <Initials>CC</Initials>
                    <Identifier Source="ORCID">http://orcid.org/0000-0002-2648-1522</Identifier>
                    <AffiliationInfo>
                        <Affiliation domain="tw" place="Taiwan">Department of Emergency Medicine, National Taiwan University Hospital, Taipei, Taiwan Department of Emergency Medicine, National Taiwan University Hospital, Yunlin Branch, Douliou, Taiwan.</Affiliation>
 (...)
```

## See also

 * Mapping NCBI/PUBMED: http://plindenbaum.blogspot.fr/2007/06/mapping-ncbipubmed.html



