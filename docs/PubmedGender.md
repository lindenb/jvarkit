# PubmedGender

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
      What kind of help
      Possible Values: [usage, markdown, xml]
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make pubmedgender
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGender.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedGender.java
)
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

## See also


 * A Simple tool to get the sex ratio in pubmed :  http://plindenbaum.blogspot.fr/2010/09/simple-tool-to-get-sex-ratio-in-pubmed.html


 

