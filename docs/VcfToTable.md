# VcfToTable

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert a vcf to a table, to ease display in the terminal


## Usage

```
Usage: vcf2table [options] Files
  Options:
    --chartsize
      google charts dimension (HTML only). Format (integer)x(integer). eg: 
      '1000x500' or (width) e.g: '1000'
    --color, --colors
      [20170808] Print Terminal ANSI colors.
      Default: false
    --format
      [20171020] output format.
      Default: text
      Possible Values: [text, html]
    --google
      use google charts (HTML only)
      Default: false
    -H, --header
      Print Header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hideAlleles
      [20170808] hide Alleles table.
      Default: false
    --hideFilters
      [20170808] hide Filters table.
      Default: false
    --hideGTypes
      [20180221] hide Genotype.Type table
      Default: false
    -g, --hideGenotypes
      Hide All genotypes
      Default: false
    -hr, --hideHomRefs
      Hide HOM_REF genotypes
      Default: false
    --hideHyperlinks
      [20191102] hide Hyperlinks table.
      Default: false
    --hideInfo
      [20170808] hide INFO table.
      Default: false
    -nc, --hideNoCalls
      Hide NO_CALL genotypes
      Default: false
    --hidePredictions
      [20170808] hide SNPEFF/VEP table.
      Default: false
    -L, -limit, --limit
      Limit the number of output variant. '-1' == ALL/No limit.
      Default: -1
    --no-html-header
      [20171023] ignore html header for HTML output.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --ped, --pedigree
      Optional Pedigree file:A pedigree is a text file delimited with tabs. No 
      header. Columns are (1) Family (2) Individual-ID (3) Father Id or '0' 
      (4) Mother Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 
      0 unaffected, 1 affected,-9 unknown  If undefined, this tool will try to 
      get the pedigree from the header.
    --url
      A custom URL for a web browser. The following words will be replaced by 
      their values: ${CHROM}, ${START}, ${END}. For example for IGV that would 
      be: 'http://localhost:60151/goto?locus=${CHROM}%3A${START}-${END}' (see 
      http://software.broadinstitute.org/software/igv/book/export/html/189) 
    --version
      print version and exit

```


## Keywords

 * vcf
 * table
 * visualization



## See also in Biostars

 * [https://www.biostars.org/p/293855](https://www.biostars.org/p/293855)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcf2table
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToTable.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToTable.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfToTableTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfToTableTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2table** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * Vcf2table : a VCF prettifier. Lindenbaum & al. 2018. figshare. [https://doi.org/10.6084/m9.figshare.5853801](https://doi.org/10.6084/m9.figshare.5853801)


## Example 

```
$ cat input.ped

FAM	M10475	0	0	1	1
FAM	M10478	0	0	2	0
FAM	M10500	M10475	M10478	2	1


$ curl -s "https://raw.githubusercontent.com/arq5x/gemini/master/test/test.region.vep.vcf" | java -jar dist/vcf2table.jar -H -p input.ped

 
INFO
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| ID              | Type    | Count | Description                                                                                                                                             |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| AC              | Integer |       | Allele count in genotypes, for each ALT allele, in the same order as listed                                                                             |
| AF              | Float   |       | Allele Frequency, for each ALT allele, in the same order as listed                                                                                      |
| AN              | Integer | 1     | Total number of alleles in called genotypes                                                                                                             |
| BaseQRankSum    | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities                                                                                       |
| CSQ             | String  |       | Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM |
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
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+

FORMAT
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
| chrX                  | 155270560 | hg19 |
| chrY                  | 59373566  | hg19 |
+-----------------------+-----------+------+

Samples
+--------+---------+--------+--------+--------+------------+
| Family | Sample  | Father | Mother | Sex    | Status     |
+--------+---------+--------+--------+--------+------------+
| FAM    | M10475  |        |        | male   | affected   |
| FAM    | M10478  |        |        | female | unaffected |
| FAM    | M10500  | M10475 | M10478 | female | affected   |
| FAM    | M128215 | M10500 |        | male   | unaffected |
+--------+---------+--------+--------+--------+------------+

>>chr1/10001/T (n 1)
 Variant
 +--------+--------------------+
 | Key    | Value              |
 +--------+--------------------+
 | CHROM  | chr1               |
 | POS    | 10001              |
 | end    | 10001              |
 | ID     | .                  |
 | REF    | T                  |
 | ALT    | TC                 |
 | QUAL   | 175.91000000000003 |
 | FILTER |                    |
 | Type   | INDEL              |
 +--------+--------------------+
 Alleles
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 | Idx | REF | Sym | Bases | Length | AC | AN | AF  | AC_affected | AC_unaffected | AC_male | AC_female |
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 | 0   | *   |     | T     | 1      | 4  | 8  | 0.5 | 2           | 1             | 1       | 2         |
 | 1   |     |     | TC    | 2      | 4  | 8  | 0.5 | 2           | 1             | 1       | 2         |
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 INFO
 +----------------+-------+----------+
 | key            | Index | Value    |
 +----------------+-------+----------+
 | AC             |       | 4        |
 | AF             |       | 0.50     |
 | AN             |       | 8        |
 | BaseQRankSum   |       | 4.975    |
 | DP             |       | 76       |
 | FS             |       | 12.516   |
 | HRun           |       | 0        |
 | HaplotypeScore |       | 218.6157 |
 | MQ             |       | 35.31    |
 | MQ0            |       | 0        |
 | MQRankSum      |       | -0.238   |
 | QD             |       | 2.31     |
 | ReadPosRankSum |       | 2.910    |
 +----------------+-------+----------+
 VEP
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 | PolyPhen                 | EXON | SIFT           | ALLELE_NUM | Gene            | SYMBOL | Protein_position | Consequence                                   | Amino_acids | Codons  | Feature         | BIOTYPE              |
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 | probably_damaging(0.956) | 8/9  | deleterious(0) | 1          | ENSG00000102967 | DHODH  | 346/395          | missense_variant                              | R/W         | Cgg/Tgg | ENST00000219240 | protein_coding       |
 |                          | 3/4  |                | 1          | ENSG00000102967 | DHODH  |                  | non_coding_exon_variant&nc_transcript_variant |             |         | ENST00000571392 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000572003 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000573843 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000573922 | processed_transcript |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  | -/193            | intron_variant                                |             |         | ENST00000574309 | protein_coding       |
 | probably_damaging(0.946) | 8/9  | deleterious(0) | 1          | ENSG00000102967 | DHODH  | 344/393          | missense_variant                              | R/W         | Cgg/Tgg | ENST00000572887 | protein_coding       |
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 Genotypes
 +---------+------+-------+----+----+-----+---------+
 | Sample  | Type | AD    | DP | GQ | GT  | PL      |
 +---------+------+-------+----+----+-----+---------+
 | M10475  | HET  | 10,2  | 15 | 10 | 0/1 | 25,0,10 |
 | M10478  | HET  | 10,4  | 16 | 5  | 0/1 | 40,0,5  |
 | M10500  | HET  | 10,10 | 21 | 7  | 0/1 | 111,0,7 |
 | M128215 | HET  | 15,5  | 24 | 0  | 0/1 | 49,0,0  |
 +---------+------+-------+----+----+-----+---------+
 TRIOS
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
 | Father-ID | Father-GT | Mother-ID | Mother-GT | Child-ID | Child-GT | Incompat. |
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
 | M10475    | 0/1       | M10478    | 0/1       | M10500   | 0/1      |           |
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
<<chr1/10001/T n 1

(...)
```

### Html output:


```
$ java -jar dist/vcf2table.jar file.vcf --color --format html > out.html
```

### Screenshots

[https://twitter.com/yokofakun/status/1067730485487366145](https://twitter.com/yokofakun/status/1067730485487366145)

![https://pbs.twimg.com/media/DtFXEhLWkAE5roc.jpg](https://pbs.twimg.com/media/DtFXEhLWkAE5roc.jpg)

[https://twitter.com/yokofakun/status/922475502933368832](https://twitter.com/yokofakun/status/922475502933368832)

![https://pbs.twimg.com/media/DM1KdWFX0AUfbxR.jpg](https://pbs.twimg.com/media/DM1KdWFX0AUfbxR.jpg)

