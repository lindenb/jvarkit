# NcbiTaxonomyToXml

Dump NCBI taxonomy tree as a hierarchical XML document


## Usage

```
Usage: ncbitaxonomy2xml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    -r, --root
      NCBI taxon root id.
      Default: 1
    --version
      print version and exit

```


## Keywords

 * taxonomy
 * ncbi
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/10327](https://www.biostars.org/p/10327)


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
$ make ncbitaxonomy2xml
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NcbiTaxonomyToXml.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NcbiTaxonomyToXml.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ncbitaxonomy2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

producing taxonomy.xml using make:

```bash
taxonomy.xml: nodes.dmp names.dmp
	java -jar dist/ncbitaxonomy2xml.jar . | xmllint --format - > $@

nodes.dmp : taxdump.tar.gz
	tar xvfz $< $@
names.dmp :taxdump.tar.gz
	tar xvfz $< $@

taxdump.tar.gz:
	curl -o $@  "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

```

output:
```xml
<?xml version="1.0" encoding="UTF-8"?>
<TaxaSet>
  <!--NcbiTaxonomyToXml by Pierre Lindenbaum PhD.-->
  <!--.-->
  <Taxon id="1">
    <TaxId>1</TaxId>
    <ScientificName>root</ScientificName>
    <TaxaSet>
      <Taxon id="131567">
        <TaxId>131567</TaxId>
        <ScientificName>cellular organisms</ScientificName>
        <ParentTaxId>1</ParentTaxId>
        <TaxaSet>
          <Taxon id="2">
            <TaxId>2</TaxId>
            <Rank>superkingdom</Rank>
            <ScientificName>Bacteria</ScientificName>
            <ParentTaxId>131567</ParentTaxId>
            <TaxaSet>
              <Taxon id="51290">
                <TaxId>51290</TaxId>
                <Rank>superphylum</Rank>
                <ScientificName>Chlamydiae/Verrucomicrobia group</ScientificName>
                <ParentTaxId>2</ParentTaxId>
                <TaxaSet>
                  <Taxon id="256845">
                    <TaxId>256845</TaxId>
                    <Rank>phylum</Rank>
                    <ScientificName>Lentisphaerae</ScientificName>
                    <ParentTaxId>51290</ParentTaxId>
                    <TaxaSet>
                      <Taxon id="278094">
                        <TaxId>278094</TaxId>
                        <ScientificName>environmental samples</ScientificName>
                        <ParentTaxId>256845</ParentTaxId>
                        <TaxaSet>
                          <Taxon id="278095">
                            <TaxId>278095</TaxId>
                            <Rank>species</Rank>
                            <ScientificName>uncultured Lentisphaerae bacterium</ScientificName>
                            <ParentTaxId>278094</ParentTaxId>
                            <TaxaSet/>
                          </Taxon>
(...)
```


