# GtfLiftOver

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

LiftOver GTF file.


## Usage

```
Usage: gtfliftover [options] Files
  Options:
  * -f, --chain
      LiftOver file.
    -x, --failed
      write  failing the liftOver here. Optional.
    -H, --header
      include gtf header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --minmatch
      lift over min-match.
      Default: 0.95
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gtf
 * liftover


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gtfliftover
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190823

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfLiftOver.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfLiftOver.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtfliftover** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```
$ gunzip -c ~/data/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz | grep ^X | head 
X	Regulatory_Build	CTCF_binding_site	100072201	100072600	.	.	.	ID=CTCF_binding_site:ENSR00000911312;bound_end=100072600;bound_start=100072201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100076001	100076400	.	.	.	ID=CTCF_binding_site:ENSR00000911313;bound_end=100076400;bound_start=100076001;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100080601	100080800	.	.	.	ID=CTCF_binding_site:ENSR00000911314;bound_end=100080800;bound_start=100080601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100099201	100099400	.	.	.	ID=CTCF_binding_site:ENSR00000911321;bound_end=100099400;bound_start=100099201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10015601	10016000	.	.	.	ID=CTCF_binding_site:ENSR00000901123;bound_end=10016000;bound_start=10015601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10018001	10018400	.	.	.	ID=CTCF_binding_site:ENSR00000901124;bound_end=10018400;bound_start=10018001;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100185401	100185800	.	.	.	ID=CTCF_binding_site:ENSR00000911326;bound_end=100185800;bound_start=100185401;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	10021601	10021800	.	.	.	ID=CTCF_binding_site:ENSR00000901125;bound_end=10021800;bound_start=10021601;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100351201	100351400	.	.	.	ID=CTCF_binding_site:ENSR00001161915;bound_end=100351400;bound_start=100351201;description=CTCF binding site;feature_type=CTCF Binding Site
X	Regulatory_Build	CTCF_binding_site	100368001	100368800	.	.	.	ID=CTCF_binding_site:ENSR00000342428;bound_end=100368800;bound_start=100368001;description=CTCF binding site;feature_type=CTCF Binding Site

$ java -jar dist/gtfliftover.jar -x jeter.boum -f ../htsjdk/src/test/resources/htsjdk/samtools/liftover/hg18ToHg19.over.chain ~/data/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz  2> /dev/null | head
chrX	Regulatory_Build	CTCF_binding_site	100185545	100185944	.	.	.	ID=CTCF_binding_site:ENSR00000911312;bound_end=100072600;bound_start=100072201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100189345	100189744	.	.	.	ID=CTCF_binding_site:ENSR00000911313;bound_end=100076400;bound_start=100076001;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100193945	100194144	.	.	.	ID=CTCF_binding_site:ENSR00000911314;bound_end=100080800;bound_start=100080601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100212545	100212744	.	.	.	ID=CTCF_binding_site:ENSR00000911321;bound_end=100099400;bound_start=100099201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10055601	10056000	.	.	.	ID=CTCF_binding_site:ENSR00000901123;bound_end=10016000;bound_start=10015601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10058001	10058400	.	.	.	ID=CTCF_binding_site:ENSR00000901124;bound_end=10018400;bound_start=10018001;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100298745	100299144	.	.	.	ID=CTCF_binding_site:ENSR00000911326;bound_end=100185800;bound_start=100185401;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	10061601	10061800	.	.	.	ID=CTCF_binding_site:ENSR00000901125;bound_end=10021800;bound_start=10021601;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100464545	100464744	.	.	.	ID=CTCF_binding_site:ENSR00001161915;bound_end=100351400;bound_start=100351201;description=CTCF binding site;feature_type=CTCF Binding Site
chrX	Regulatory_Build	CTCF_binding_site	100481345	100482144	.	.	.	ID=CTCF_binding_site:ENSR00000342428;bound_end=100368800;bound_start=100368001;description=CTCF binding site;feature_type=CTCF Binding Site
```

when the GTF contains a gene, all sub-section MUST be re-mapped.

```
$ gunzip -c ~/src/jvarkit-git/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | head
22	ensembl_havana	exon	41697526	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00001307363"; exon_version "4"; tag "basic";
22	ensembl_havana	five_prime_utr	41697526	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	ensembl_havana	gene	41697526	41756151	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
22	ensembl_havana	transcript	41697526	41756151	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	havana	exon	41697719	41697776	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00001942555"; exon_version "1";
22	havana	transcript	41697719	41732847	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1";
22	ensembl_havana	exon	41716659	41716717	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "2"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00003628979"; exon_version "1"; tag "basic";
22	ensembl_havana	five_prime_utr	41716659	41716664	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
22	havana	exon	41716659	41716717	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "2"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00003530265"; exon_version "1";
22	ensembl	CDS	41716665	41716717	.	+	0	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000351589"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; protein_id "ENSP00000263243"; protein_version "5"; tag "basic";

$ java -jar dist/gtfliftover.jar -f ../htsjdk/src/test/resources/htsjdk/samtools/liftover/hg18ToHg19.over.chain ~/src/jvarkit-git/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz  --header 2> /dev/null  | head
chr22	ensembl_havana	exon	43367582	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; exon_id "ENSE00001307363"; exon_version "4"; tag "basic";
chr22	ensembl_havana	five_prime_utr	43367582	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
chr22	ensembl_havana	gene	43367582	43426207	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
chr22	ensembl_havana	transcript	43367582	43426207	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000352645"; transcript_version "4"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14013"; havana_transcript "OTTHUMT00000320696"; havana_transcript_version "1"; tag "basic";
chr22	havana	exon	43367775	43367832	.	+	.	gene_id "ENSG00000100403"; gene_version "10"; transcript_id "ENST00000486331"; transcript_version "1"; exon_number "1"; gene_name "ZC3H7B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ZC3H7B-002"; transcript_source "havana"; transcript_biotype "retained_intron"; havana_transcript "OTTHUMT00000320697"; havana_transcript_version "1"; exon_id "ENSE00001942555"; exon_version "1";
```

