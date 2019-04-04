# NcbiGeneDump

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Dump XML results from a list of gene using NCBI/Eutils


## Usage

```
Usage: ncbigenedump [options] Files
  Options:
    --abort
      Abort with error if a gene was not found with ncbi-esearch.
      Default: false
    -C, --custom
      Custom annotation file. See the main documentation.
    -e, --email
      optional user email
    -L, -G, --list, --genes
      File containing a list of genes, can be a gene name or a ncbi gene id, 
      one per line.
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
    --seconds
      wait 'n' seconds between each calls.
      Default: 2
    -skip, --skip
      Optional set of elements names to be ignored in the output. Spaces or 
      comma separated. .eg: 'Gene-track '
      Default: <empty string>
    --stdin
      read list of genes from stdin.
      Default: false
    -T, --taxon
      taxon id.
      Default: 9606
    --version
      print version and exit

```


## Keywords

 * ncbi
 * gene
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew ncbigenedump
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NcbiGeneDump.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/NcbiGeneDump.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/NcbiGeneDumpTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/NcbiGeneDumpTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ncbigenedump** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Annotation files

a custom annotation file can be a xml or a plain text file.

Annotation will be inserted in a `<custom-annotation>` tag (not in the official dtd) under `<Entrezgene>`

### XML files

The program will insert the first XML tag containing a XML attribute `ncbi-gene-id` corresponding to the current NCBI gene identifier.

### Plain text files

block of lines starting after `## NCBI Gene ID.` followed by the current NCBI gene identifer will be inserted.



## Examples

### Example 1

search two genes and convert to markdown using the following XSLT stylesheet.

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'  version='1.0' >
<xsl:output method="text" encoding="UTF-8"/>

<xsl:template match="Entrezgene-Set">
<xsl:apply-templates select="Entrezgene"/>
</xsl:template>

<xsl:template match="Entrezgene">
<xsl:variable name="geneId" select="Entrezgene_track-info/Gene-track/Gene-track_geneid"/>
<xsl:variable name="locus" select="Entrezgene_gene/Gene-ref/Gene-ref_locus"/>
<xsl:text>**</xsl:text>
<xsl:value-of select="$locus"/>
<xsl:text>** : </xsl:text>
<xsl:value-of select="Entrezgene_summary"/>

<xsl:if test="custom-annotation">
<xsl:text>
CUSTOM-ANNOTATIONS: </xsl:text><xsl:value-of select="normalize-space(custom-annotation)"/>
</xsl:if>
<xsl:text>

</xsl:text>
</xsl:template>

</xsl:stylesheet>
```


```
$ java -jar dist/ncbigenedump.jar SCN5A NOTCH2 |\
 	xsltproc stylesheet.xsl  -
```

**NOTCH2** : This gene encodes a member of the Notch family. Members of this Type 1 transmembrane protein family share structural characteristics including an extracellular domain consisting of multiple epidermal growth factor-like (EGF) repeats, and an intracellular domain consisting of multiple, different domain types. Notch family members play a role in a variety of developmental processes by controlling cell fate decisions. The Notch signaling network is an evolutionarily conserved intercellular signaling pathway which regulates interactions between physically adjacent cells. In Drosophilia, notch interaction with its cell-bound ligands (delta, serrate) establishes an intercellular signaling pathway that plays a key role in development. Homologues of the notch-ligands have also been identified in human, but precise interactions between these ligands and the human notch homologues remain to be determined. This protein is cleaved in the trans-Golgi network, and presented on the cell surface as a heterodimer. This protein functions as a receptor for membrane bound ligands, and may play a role in vascular, renal and hepatic development. Two transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jan 2011]

**SCN5A** : The protein encoded by this gene is an integral membrane protein and tetrodotoxin-resistant voltage-gated sodium channel subunit. This protein is found primarily in cardiac muscle and is responsible for the initial upstroke of the action potential in an electrocardiogram. Defects in this gene are a cause of long QT syndrome type 3 (LQT3), an autosomal dominant cardiac disease. Alternative splicing results in several transcript variants encoding different isoforms. [provided by RefSeq, Jul 2008]

### Example

The xml annotation file:

```html
<html xmlns="http://www.w3.org/1999/xhtml"><head><title>My custom annots</title></head><body>
<div ncbi-gene-id="6331"><h2>SCN5A</h2><p>This is the gene, Matilde.</p></div>
<div ncbi-gene-id="4853"><h2>NOTCH2</h2><p>Hajdu-Cheney syndrome</p></div>
</body></html>
```

```
java -jar dist/ncbigenedump.jar -C annot.html NOTCH2 SCN5A | xmllint --format -

(...)
      </Xtra-Terms>
    </Entrezgene_xtra-properties>
    <custom-annotation>
      <div ncbi-gene-id="6331">
        <h2>SCN5A</h2>
        <p>This is the gene, Matilde.</p>
      </div>
    </custom-annotation>
  </Entrezgene>
</Entrezgene-Set>

```

### Example

The text annotation file `annot.md`:

```markdown
# My Annotations

last updated 2018-08-29

## NCBI Gene ID. 4853 NOTCH2

encodes a member of the Notch family

## NCBI Gene ID. 6331 SCN5A

See also SCN10A
```

```
java -jar dist/ncbigenedump.jar -C annot.md NOTCH2 SCN5A | xmllint --format -

(...)
        <Xtra-Terms_tag>PROP</Xtra-Terms_tag>
        <Xtra-Terms_value>phenotype</Xtra-Terms_value>
      </Xtra-Terms>
    </Entrezgene_xtra-properties>
    <custom-annotation>
encodes a member of the Notch family

</custom-annotation>
  </Entrezgene>
</Entrezgene-Set>
```

### Example

custom annotation and XSLT:

```
$ java -jar dist/ncbigenedump.jar -C annot.md NOTCH2 SCN5A | xsltproc tansform.xsl -
```

output:

**NOTCH2** : This gene encodes a member of the Notch family. Members of this Type 1 transmembrane protein family share structural characteristics including an extracellular domain consisting of multiple epidermal growth factor-like (EGF) repeats, and an intracellular domain consisting of multiple, different domain types. Notch family members play a role in a variety of developmental processes by controlling cell fate decisions. The Notch signaling network is an evolutionarily conserved intercellular signaling pathway which regulates interactions between physically adjacent cells. In Drosophilia, notch interaction with its cell-bound ligands (delta, serrate) establishes an intercellular signaling pathway that plays a key role in development. Homologues of the notch-ligands have also been identified in human, but precise interactions between these ligands and the human notch homologues remain to be determined. This protein is cleaved in the trans-Golgi network, and presented on the cell surface as a heterodimer. This protein functions as a receptor for membrane bound ligands, and may play a role in vascular, renal and hepatic development. Two transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jan 2011]

CUSTOM-ANNOTATIONS: encodes a member of the Notch family

**SCN5A** : The protein encoded by this gene is an integral membrane protein and tetrodotoxin-resistant voltage-gated sodium channel subunit. This protein is found primarily in cardiac muscle and is responsible for the initial upstroke of the action potential in an electrocardiogram. Defects in this gene are a cause of long QT syndrome type 3 (LQT3), an autosomal dominant cardiac disease. Alternative splicing results in several transcript variants encoding different isoforms. [provided by RefSeq, Jul 2008]

CUSTOM-ANNOTATIONS: See also SCN10A


## Screenshots


https://twitter.com/yokofakun/status/1034825482808819713

![https://twitter.com/yokofakun/status/1034825482808819713](https://pbs.twimg.com/media/DlxwR3fXoAIFKO2.jpg)

