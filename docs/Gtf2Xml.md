# Gtf2Xml

Convert GTF to XML


## Usage

```
Usage: gtf2xml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tmpDir
      mp directory
      Default: /tmp
    --version
      print version and exit

```


## Keywords

 * xml
 * gtf


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
$ make gtf2xml
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/Gtf2Xml.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/Gtf2Xml.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtf2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$ curl  "ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz" | gunzip -c |\
 head -n 10 | java -jar dist/gtf2xml.jar | xmllint --format -

```
output: 
```xml
<?xml version="1.0" encoding="UTF-8"?>
<gtf gff-version="3" genome-build="GRCh38.p3" genome-version="GRCh38" genome-date="2013-12" genome-build-accession="NCBI:GCA_000001405.18" genebuild-last-updated="2015-06">
  <feature chrom="1" start="11869" end="14409" strand="+" source="havana" type="gene">
    <attributes>
      <ID>gene:ENSG00000223972</ID>
      <Name>DDX11L1</Name>
      <biotype>transcribed_unprocessed_pseudogene</biotype>
      <description>DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1 [Source:HGNC Symbol;Acc:HGNC:37102]</description>
      <gene_id>ENSG00000223972</gene_id>
      <havana_gene>OTTHUMG00000000961</havana_gene>
      <havana_version>2</havana_version>
      <logic_name>havana</logic_name>
      <version>5</version>
    </attributes>
  </feature>
  <feature chrom="1" start="14404" end="29570" strand="-" source="havana" type="gene">
    <attributes>
      <ID>gene:ENSG00000227232</ID>
      <Name>WASH7P</Name>
      <biotype>unprocessed_pseudogene</biotype>
      <description>WAS protein family homolog 7 pseudogene [Source:HGNC Symbol;Acc:HGNC:38034]</description>
      <gene_id>ENSG00000227232</gene_id>
      <havana_gene>OTTHUMG00000000958</havana_gene>
      <havana_version>1</havana_version>
      <logic_name>havana</logic_name>
      <version>5</version>
    </attributes>
  </feature>
  <feature chrom="1" start="17369" end="17436" strand="-" source="ensembl" type="miRNA_gene">
    <attributes>
      <ID>gene:ENSG00000278267</ID>
      <Name>MIR6859-1</Name>
      <biotype>miRNA</biotype>
      <description>microRNA 6859-1 [Source:HGNC Symbol;Acc:HGNC:50039]</description>
      <gene_id>ENSG00000278267</gene_id>
      <logic_name>ncrna</logic_name>
      <version>1</version>
    </attributes>
  </feature>
  <feature chrom="1" start="29554" end="31109" strand="+" source="havana" type="lincRNA_gene">
    <attributes>
      <ID>gene:ENSG00000243485</ID>
      <Name>RP11-34P13.3</Name>
      <biotype>lincRNA</biotype>
      <gene_id>ENSG00000243485</gene_id>
      <havana_gene>OTTHUMG00000000959</havana_gene>
      <havana_version>2</havana_version>
      <logic_name>havana</logic_name>
      <version>3</version>
    </attributes>
  </feature>
  <attributes>
    <attribute>havana_gene</attribute>
    <attribute>Name</attribute>
    <attribute>havana_version</attribute>
    <attribute>logic_name</attribute>
    <attribute>description</attribute>
    <attribute>biotype</attribute>
    <attribute>ID</attribute>
    <attribute>gene_id</attribute>
    <attribute>version</attribute>
  </attributes>
  <types>
    <type>miRNA_gene</type>
    <type>gene</type>
    <type>lincRNA_gene</type>
  </types>
  <sources>
    <source>ensembl</source>
    <source>havana</source>
  </sources>
  <dict>
    <chrom name="1" length="31109"/>
  </dict>
</gtf>
```

 

