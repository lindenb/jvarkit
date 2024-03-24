# NcbiTaxonomyToXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Dump NCBI taxonomy tree as a hierarchical XML document or as a table


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar ncbitaxonomy2xml  [options] Files

Usage: ncbitaxonomy2xml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -r, --root
      NCBI taxon root id.
      Default: 1
    -t, --tabular
      produce tabular output instead of xml
      Default: false
    --version
      print version and exit

```


## Keywords

 * taxonomy
 * ncbi
 * xml



## See also in Biostars

 * [https://www.biostars.org/p/10327](https://www.biostars.org/p/10327)



## Creation Date

20120320

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/taxonomy/NcbiTaxonomyToXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/taxonomy/NcbiTaxonomyToXml.java)


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


