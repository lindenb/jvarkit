# VcfBiomart

BiomartQueries with VCF


## Usage

```
Usage: vcfbiomart [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -C
       (int) column index (1-based) for chromosome . Optional
      Default: -1
    -E
       (int) column index (1-based) for end . Optional
      Default: -1
    -S
       (int) column index (1-based) for start . Optional
      Default: -1
    -T
       (string) VCF output tag.
      Default: BIOMART
    -X
       (XML-file) XML biomart template.
    -n
       (int) batch size.);
      Default: 200
    -u
       (url) biomart service url
      Default: http://www.biomart.org/biomart/martservice/result

```


## Keywords

 * vcf
 * ensembl
 * biomart


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
$ make vcfbiomart
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbiomart/VcfBiomart.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbiomart/VcfBiomart.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbiomart** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

Let's annotate a VCF with the gene-ncbi/gene-start/gene-end/gene-ncbiGeneId .

From Ensembl, we've downloaded the following XML file:

```xml
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```
This file is used  as follow:

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  -C 1 -S 2 -E 3 |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:chromosome_name|start_position|end_position|entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=1|145209119|145291972|388677,1|145209145|145319796|	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=10|1095478|1178237|22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=10|135204338|135233999|,10|135207598|135234811|92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
###XML attributes: @chrom/@start/@end.
By adding the attributes chrom=true, start=true end=true the column for chrom/start/end can be specified in the XML
```XML
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true"/>
		<Attribute name = "start_position" start="true"/>
		<Attribute name = "end_position" end="true" />
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```
```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:chromosome_name|start_position|end_position|entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=1|145209119|145291972|388677,1|145209145|145319796|	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=10|1095478|1178237|22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=10|135204338|135233999|,10|135207598|135234811|92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
###Hiding columns: @visible=false
By adding the attribute visible=false, some columns can be removed from the result.
```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true" visible="false"/>
		<Attribute name = "start_position" start="true" visible="false"/>
		<Attribute name = "end_position" end="true" visible="false"/>
		<Attribute name = "entrezgene" />
	</Dataset>
</Query>
```

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=.,Type=String,Description="Biomart query. Format:entrezgene">
1	145273345	.	T	C	289.85	.	NCBIGENE=388677	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE=22884	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE=92170	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```
if all the attributes are set to visible=false, the INFO is set to 'Flag'
```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "chromosome_name" chrom="true" visible="false"/>
		<Attribute name = "start_position" start="true" visible="false"/>
		<Attribute name = "end_position" end="true" visible="false"/>
		<Attribute name = "entrezgene" visible="false"/>
	</Dataset>
</Query>
```

```bash
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  sed 's/chr//g' |\
 java -jar dist/vcfbiomart.jar -X mart.xml -T NCBIGENE  |\
 grep NCBIGENE

##INFO=<ID=NCBIGENE,Number=0,Type=Flag,Description="Biomart query.">
1	145273345	.	T	C	289.85	.	NCBIGENE	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
10	1142208	.	T	C	3404.30	.	NCBIGENE	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
10	135210791	.	T	C	65.41	.	NCBIGENE	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
```

