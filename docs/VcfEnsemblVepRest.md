# VcfEnsemblVepRest

Annotate a VCF with ensembl REST API


## Usage

```
Usage: vcfensemblvep [options] Files
  Options:
    -x, --base64
      save whole XML document as xml base 64
      Default: false
    -n, --batchSize
      batch size
      Default: 100
    -e, --extension
      Path extension
      Default: /vep/homo_sapiens/region
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -s, --server
      REST server
      Default: http://grch37.rest.ensembl.org
    -T, --tee
      'Tee' xml response to stderr
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * rest
 * ensembl


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
$ make vcfensemblvep
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ensembl/VcfEnsemblVepRest.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ensembl/VcfEnsemblVepRest.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfensemblvep** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




## Example

```bash
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c  | java -jar dist/vcfensemblvep.jar | grep -v '^#' | cut -f 1,2,4,5,8


1	10583	G	A	VEPTRCSQ=processed_transcript|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|A|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|A|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|A|SO:0001631,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|A|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|A|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|A|SO:0001632
1	10611	C	G	VEPTRCSQ=processed_transcript|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|G|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|G|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|G|SO:0001631,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|G|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|G|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|G|SO:0001632
1	13302	C	T	VEPTRCSQ=processed_transcript|550|550|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|T|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|T|SO:0001632,transcribed_unprocessed_pseudogene|342|342|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|T|SO:0001792&SO:0001619,transcribed_unprocessed_pseudogene|543|543|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|T|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|T|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|T|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|T|SO:0001632
1	13327	G	C	VEPTRCSQ=processed_transcript|575|575|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|C|SO:0001632,transcribed_unprocessed_pseudogene|367|367|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|C|SO:0001792&SO:0001619,transcribed_unprocessed_pseudogene|568|568|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|C|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|C|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|C|SO:0001632
1	13957	TC	T	VEPTRCSQ=processed_transcript|1206|1206|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147||SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675||SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305||SO:0001632,transcribed_unprocessed_pseudogene|1199|1199|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476||SO:0001632,transcribed_unprocessed_pseudogene|1032|1032|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504||SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562||SO:0001632
1	13980	T	C	VEPTRCSQ=processed_transcript|1228|1228|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|C|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|C|SO:0001632,transcribed_unprocessed_pseudogene|1221|1221|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|C|SO:0001632,transcribed_unprocessed_pseudogene|1054|1054|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|C|SO:0001632
1	30923	G	T	VEPTRCSQ=lincRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000473358|T|SO:0001627&SO:0001619,lincRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000469289|T|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|T|SO:0001631,lincRNA|||||ENSG00000237613|FAM138A|HGNC|32334|-1|ENST00000417324|T|SO:0001632,miRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000607096|T|SO:0001632,lincRNA|||||ENSG00000237613|FAM138A|HGNC|32334|-1|ENST00000461467|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|T|SO:0001631
1	46402	C	CTGT	.
1	47190	G	GA	.
1	51476	T	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	51479	T	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001631
1	51914	T	G	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|G|SO:0001631
1	51935	C	T	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|T|SO:0001631
1	51954	G	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	52058	G	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	52144	T	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001631
1	52185	TTAA	T	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647||SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857||SO:0001631
1	52238	T	G	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|G|SO:0001631
[INFO/VcfEnsemblVepRest] 2015-03-19 17:19:04 "done: N=20"
[INFO/VcfEnsemblVepRest] 2015-03-19 17:19:04 "Number of Variants:0"
1	53234	CAT	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647||SO:0001627&SO:0001619,unprocessed_pseudogene|763|764|||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857||SO:0001792&SO:0001619
1	54353	C	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001632

```


