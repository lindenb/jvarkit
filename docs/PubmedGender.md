# PubmedGender

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Add gender-related attributes in the Author tag of pubmed xml. 


## Usage

```
Usage: pubmedgender [options] Files
  Options:
  * -d, --database
      REQUIRED: A comma delimited file containing the following columns: 1) 
      Name 2) sex (M/F) 3) Score. See 
      http://cpansearch.perl.org/src/EDALY/Text-GenderFromName-0.33/GenderFromName.pm 
      or https://www.ssa.gov/oact/babynames/names.zip
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
 * gender
 * ncbi
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pubmedgender
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGender.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGender.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGenderTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGenderTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedgender** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Building the database

```
$ wget -O jeter.zip "https://www.ssa.gov/oact/babynames/names.zip"
$ unzip -t jeter.zip | tail
    testing: yob2009.txt              OK
    testing: yob2010.txt              OK
    testing: yob2011.txt              OK
    testing: yob2012.txt              OK
    testing: yob2013.txt              OK
    testing: yob2014.txt              OK
    testing: yob2015.txt              OK
    testing: yob1880.txt              OK
    testing: NationalReadMe.pdf       OK
No errors detected in compressed data of jeter.zip.
$ unzip -p jeter.zip yob2015.txt &gt; database.csv
```

## Example

```
$ java -jar dist/pubmeddump.jar "Lindenbaum[Author] Nantes" 2> /dev/null  | java -jar dist/pubmedgender.jar  -d jeter.csv 2> /dev/null | grep Lindenbaum -A 2 -B 1
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
```

## Example


Get an histogram for gender ratio in pubmed:

```
$ java -jar dist/pubmeddump.jar  '("bioinformatics"[Journal]) AND ("2018-01-01"[Date - Publication] : "2018-07-06"[Date - Publication]) ' |\
	java -jar dist/pubmedgender.jar -d ./yob2017.txt  |\
	java -jar dist/xsltstream.jar -t ~/tr.xsl -n "PubmedArticle" |\
	tr ";" "\n" | sort | uniq -c |\
	java -jar dist/simpleplot.jar -su -t SIMPLE_HISTOGRAM  --title "Authors in Bioinformatics 2018"
```

with tr.xsl:

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet 
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:j='http://www.ibm.com/xmlns/prod/2009/jsonx'
	version='1.0'>
<xsl:output method="text"/>


<xsl:template match="/">
<xsl:apply-templates select="PubmedArticle"/>
</xsl:template>

<xsl:template match="PubmedArticle">
<xsl:apply-templates select="MedlineCitation/Article/AuthorList/Author[@female or @male]"/>
</xsl:template>

<xsl:template match="Author">
<xsl:choose>
	<xsl:when test="@female and @male and number(@female) &gt; number(@male) ">
		<xsl:text>FEMALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@female and @male and number(@female) &lt; number(@male) ">
		<xsl:text>MALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@female">
		<xsl:text>FEMALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@male">
		<xsl:text>MALE;</xsl:text>
	</xsl:when>
</xsl:choose>
</xsl:template>

</xsl:stylesheet>
```


## See also


 * A Simple tool to get the sex ratio in pubmed :  http://plindenbaum.blogspot.fr/2010/09/simple-tool-to-get-sex-ratio-in-pubmed.html


 
