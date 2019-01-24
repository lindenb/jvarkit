# VcfBiomart

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

BiomartQueries with VCF


## Usage

```
Usage: vcfbiomart [options] Files
  Options:
    -C, --contig, --chrom
      @name attribute <Filter> in the <Query> xml document.
      Default: chromosome_name
    -E, --end
       (int) column index (1-based) for end . Optional
      Default: end
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -label, --label, --labels
      Add the field label in the INFO attribute 'label1|value1|label2|value2'
      Default: false
    -u, --url, --mart
       (url) biomart service url. See 
      http://grch37.ensembl.org/info/data/biomart/biomart_restful.html 
      Default: http://grch37.ensembl.org/biomart/martservice
    -o, --output
      Output file. Optional . Default: stdout
    -S, --start
       (int) column index (1-based) for start . Optional
      Default: start
    -T, --tag
       (string) VCF output tag.
      Default: BIOMART
    -tee, --tee
      'Tee' response to stderr
      Default: false
    --version
      print version and exit
  * -X, --xml
       (XML-file) XML biomart template. <Query>...</Query>

```


## Keywords

 * vcf
 * ensembl
 * biomart
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfbiomart
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ensembl/VcfBiomart.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ensembl/VcfBiomart.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/ensembl/VcfBiomartTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/ensembl/VcfBiomartTest.java)


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


## History

* rewritten 2018-02-19

## Example

the XML query:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "end" value = "10000000"/>
		<Filter name = "start" value = "1000000"/>
		<Filter name = "chromosome_name" value = "4"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "external_transcript_name" />
		<Attribute name = "transcription_start_site" />
		<Attribute name = "transcript_start" />
		<Attribute name = "transcript_end" />
	</Dataset>
</Query>
```

the XML file above contains three fields that will be used/replaced to query the position of a variant:

```xml
	<Filter name = "end" value = "10000000"/>
	<Filter name = "start" value = "1000000"/>
	<Filter name = "chromosome_name" value = "4"/>
```

running the query:

```xml
java  dist/vcfbiomart.jar -X query.xml input.vcf
```

```
##INFO=<ID=BIOMART,Number=.,Type=String,Description="Biomart query. Format: ensembl_gene_id|ensembl_transcript_id|start_position|end_position|external_transcript_name|transcription_start_site|transcript_start|transcript_end">
(...)
(...)A|||||4113|;BIOMART=ENSG00000183873|ENST00000414099|38589548|38691164|SCN5A-004|38674840|38589548|38674840,ENSG00000183873|ENST00000423572|38589548|38691164|SCN5A-003|38674853|38589553|38674853,ENSG00000183873|ENST00000413689|38589548|38691164|SCN5A-001|38691163|38589553|38691163,ENSG00000183873|ENST00000333535|38589548|38691164|SCN5A-014|38691119|38589557|38691119,ENSG00000183873|ENST00000455624|38589548|38691164|SCN5A-002|38674823|38590619|38674823,ENSG00000183873|ENST00000450102|38589548|38691164|SCN5A-008|38674840|38591459|38674840,ENSG00000183873|ENST00000449557|38589548|38691164|SCN5A-010|38674807|38591812|38674807,ENSG00000183873|ENST00000464652|38589548|38691164|SCN5A-011|38596040|38594472|38596040,ENSG00000183873|ENST00000491944|38589548|38691164|SCN5A-005|38691164|38655519|38691164,ENSG00000183873|ENST00000476683|38589548|38691164|SCN5A-013|38674711|38671333|38674711,ENSG00000183873|ENST00000327956|38589548|38691164|SCN5A-006|38687267|38674602|38687267,ENSG00000183873|ENST00000451551|38589548|38691164|SCN5A-203|38691164|38589553|38691164,ENSG00000183873|ENST00000443581|38589548|38691164|SCN5A-202|38691163|38589553|38691163,ENSG00000183873|ENST00000425664|38589548|38691164|SCN5A-201|38691163|38589553|38691163;(...)
(...)
```

