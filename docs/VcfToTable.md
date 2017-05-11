# VcfToTable


## Usage

```
Usage: vcf2table [options] Files
  Options:
    -h, --help
      print help and exits
    -g, --hideGenotypes
      Hide All genotypes
      Default: false
    -hr, --hideHomRefs
      Hide HOM_REF genotypes
      Default: false
    -nc, --hideNoCalls
      Hide NO_CALL genotypes
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exits
    -H
      Print Header
      Default: false

```


## Description

convert a vcf to a table, to ease display in the terminal


## Keywords

 * vcf
 * table
 * visualization


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
$ make vcf2table
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

https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToTable.java

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2table** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> http://dx.doi.org/10.6084/m9.figshare.1425030


## Example 

```
$ curl -s "https://raw.githubusercontent.com/arq5x/gemini/master/test/test.region.vep.vcf" | java -jar dist/vcf2table.jar -H

INFO
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| ID              | Type    | Count | Description                                                                                                                                             |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| AC              | Integer |       | Allele count in genotypes, for each ALT allele, in the same order as listed                                                                             |
| AF              | Float   |       | Allele Frequency, for each ALT allele, in the same order as listed                                                                                      |
| AN              | Integer | 1     | Total number of alleles in called genotypes                                                                                                             |
| BaseQRankSum    | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities                                                                                       |
| DP              | Integer | 1     | Approximate read depth; some reads may have been filtered                                                                                               |
| DS              | Flag    | 0     | Were any of the samples downsampled?                                                                                                                    |
| Dels            | Float   | 1     | Fraction of Reads Containing Spanning Deletions                                                                                                         |
| FS              | Float   | 1     | Phred-scaled p-value using Fisher's exact test to detect strand bias                                                                                    |
| HRun            | Integer | 1     | Largest Contiguous Homopolymer Run of Variant Allele In Either Direction                                                                                |
| HaplotypeScore  | Float   | 1     | Consistency of the site with at most two segregating haplotypes                                                                                         |
| InbreedingCoeff | Float   | 1     | Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation                       |
| MQ              | Float   | 1     | RMS Mapping Quality                                                                                                                                     |
| MQ0             | Integer | 1     | Total Mapping Quality Zero Reads                                                                                                                        |
| MQRankSum       | Float   | 1     | Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities                                                                               |
| QD              | Float   | 1     | Variant Confidence/Quality by Depth                                                                                                                     |
| ReadPosRankSum  | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias                                                                                   |
| CSQ             | String  |       | Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+

INFO
+----+---------+-------+----------------------------------------------------------------------------------------+
| ID | Type    | Count | Description                                                                            |
+----+---------+-------+----------------------------------------------------------------------------------------+
| AD | Integer |       | Allelic depths for the ref and alt alleles in the order listed                         |
| DP | Integer | 1     | Approximate read depth (reads with MQ=255 or with bad mates are filtered)              |
| GQ | Integer | 1     | Genotype Quality                                                                       |
| GT | String  | 1     | Genotype                                                                               |
| PL | Integer |       | Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification |
+----+---------+-------+----------------------------------------------------------------------------------------+


Dict
+-----------------------+-----------+------+
| Name                  | Length    | AS   |
+-----------------------+-----------+------+
| chr1                  | 249250621 | hg19 |
(...)
| chr10                 | 135534747 | hg19 |
| chrY                  | 59373566  | hg19 |
+-----------------------+-----------+------+

>>chr1/10001/T 1
Variant
+-------+--------------------+
| Key   | Value              |
+-------+--------------------+
| CHROM | chr1               |
| POS   | 10001              |
| end   | 10001              |
| ID    | .                  |
| REF   | T                  |
| ALT   | TC                 |
| QUAL  | 175.91000000000003 |
+-------+--------------------+
Alleles
+-----+-----+-----+-------+--------+
| Idx | REF | Sym | Bases | Length |
+-----+-----+-----+-------+--------+
| 0   | *   |     | T     | 1      |
| 1   |     |     | TC    | 2      |
+-----+-----+-----+-------+--------+
INFO
+----------------+-------+----------------------------------------------------------------------------------------------------------+
| key            | Index | Value                                                                                                    |
+----------------+-------+----------------------------------------------------------------------------------------------------------+
| AC             |       | 4                                                                                                        |
| AF             |       | 0.50                                                                                                     |
| AN             |       | 8                                                                                                        |
| BaseQRankSum   |       | 4.975                                                                                                    |
| CSQ            | 1     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000456328|||||processed_transcript|1               |
| CSQ            | 2     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000488147|||||unprocessed_pseudogene|1            |
| CSQ            | 3     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000541675|||||unprocessed_pseudogene|1            |
| CSQ            | 4     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000450305|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 5     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000515242|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 6     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000538476|||||unprocessed_pseudogene|1            |
| CSQ            | 7     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000518655|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 8     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000438504|||||unprocessed_pseudogene|1            |
| CSQ            | 9     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000423562|||||unprocessed_pseudogene|1            |
| DP             |       | 76                                                                                                       |
| FS             |       | 12.516                                                                                                   |
| HRun           |       | 0                                                                                                        |
| HaplotypeScore |       | 218.6157                                                                                                 |
| MQ             |       | 35.31                                                                                                    |
| MQ0            |       | 0                                                                                                        |
| MQRankSum      |       | -0.238                                                                                                   |
| QD             |       | 2.31                                                                                                     |
| ReadPosRankSum |       | 2.910                                                                                                    |
+----------------+-------+----------------------------------------------------------------------------------------------------------+
VEP
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
| PolyPhen | EXON | SIFT | ALLELE_NUM | Gene            | SYMBOL  | Protein_position | Consequence             | Amino_acids | Codons | Feature         | BIOTYPE                            |
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000456328 | processed_transcript               |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000488147 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000541675 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000450305 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000515242 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000538476 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000518655 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000438504 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000423562 | unprocessed_pseudogene             |
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
Genotypes
+---------+------+-------+----+-------+-----+---------+
| Sample  | Type | AD    | DP | GQ    | GT  | PL      |
+---------+------+-------+----+-------+-----+---------+
| M10475  | HET  | 10,2  | 15 | 10.41 | 0/1 | 25,0,10 |
| M10478  | HET  | 10,4  | 16 | 5.30  | 0/1 | 40,0,5  |
| M10500  | HET  | 10,10 | 21 | 7.48  | 0/1 | 111,0,7 |
| M128215 | HET  | 15,5  | 24 | 0.26  | 0/1 | 49,0,0  |
+---------+------+-------+----+-------+-----+---------+
<<1
(...)
```

