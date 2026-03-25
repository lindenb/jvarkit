# GwasCatalogToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert gwas catalog TSV association file to BED using a gtf file.


## Usage

```
Usage: java -jar dist/gwascatalog2bed.jar  [options] Files
Usage: gwascatalog2bed [options] Files
  Options:
  * -G, --gtf
      GTF file
    --header
      print header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --pos, --position
      use CHR_ID/CHR_POS to get a position. Beware, check the build version.
      Default: false
    --version
      print version and exit

```


## Keywords

 * bed
 * gwas
 * gwascatalog


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gwascatalog2bed
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20260325

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gwascat2bed/GwasCatalogToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gwascat2bed/GwasCatalogToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gwascatalog2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
 $ unzip -p  "gwas-catalog.zip "gwas-catalog-download-associations-alt-full.tsv" |\
 	java -jar dist/jvarkit.jar gwascatalog2bed --gtf  hs37d5.gtf.gz --header 2> /dev/null  |\
 	verticalize 

>>> 2
$1                      #contig : chr18
$2                        start : 11981023
$3                          end : 12030876
$4        DATE ADDED TO CATALOG : 2008-06-16
$5                     PUBMEDID : 17434096
$6                 FIRST AUTHOR : Matarin M
$7                         DATE : 2007-05-06
$8                      JOURNAL : Lancet Neurol
$9                         LINK : www.ncbi.nlm.nih.gov/pubmed/17434096
$10                       STUDY : A genome-wide genotyping study in patients with ischaemic stroke: initial analysis and data release.
$11               DISEASE/TRAIT : Stroke
$12         INITIAL SAMPLE SIZE : 249 European ancestry cases, 268 European ancestry controls
$13     REPLICATION SAMPLE SIZE : NA
$14                      REGION : 18p11.21
$15                      CHR_ID : 18
$16                     CHR_POS : 11987273
$17            REPORTED GENE(S) : IMPA2
$18                 MAPPED_GENE : IMPA2
$19            UPSTREAM_GENE_ID : 
$20          DOWNSTREAM_GENE_ID : 
$21                SNP_GENE_IDS : ENSG00000141401
$22      UPSTREAM_GENE_DISTANCE : 
$23    DOWNSTREAM_GENE_DISTANCE : 
$24   STRONGEST SNP-RISK ALLELE : rs7506045-?
$25                        SNPS : rs7506045
$26                      MERGED : 0
$27              SNP_ID_CURRENT : 7506045
$28                     CONTEXT : intron_variant
$29                  INTERGENIC : 0
$30       RISK ALLELE FREQUENCY : 0.10
$31                     P-VALUE : 7E-7
$32                 PVALUE_MLOG : 6.154901959985743
$33              P-VALUE (TEXT) : 
$34                  OR or BETA : 5.39
$35               95% CI (TEXT) : [2.77-10.5]
$36  PLATFORM [SNPS PASSING QC] : Illumina [408803]
$37                         CNV : N
$38                MAPPED_TRAIT : stroke
$39            MAPPED_TRAIT_URI : http://www.ebi.ac.uk/efo/EFO_0000712
$40             STUDY ACCESSION : GCST000032
$41       GENOTYPING TECHNOLOGY : Genome-wide genotyping array
<<< 2
```


