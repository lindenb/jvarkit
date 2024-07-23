# VCFPolyX

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Number of repeated REF bases around POS.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfpolyx  [options] Files

Usage: vcfpolyx [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -n, --filter, --max-repeats
      if number of repeated bases is greater or equal to 'n' set a FILTER = 
      (tag) 
      Default: -1
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --skip-filtered
      Don't spend some time to calculate the tag if the variant is FILTERed
      Default: false
    -t, --tag
      Tag used in INFO and FILTER columns.
      Default: POLYX
    --version
      print version and exit

```


## Keywords

 * vcf
 * repeat



## Creation Date

20200930

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfpolyx/VCFPolyX.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfpolyx/VCFPolyX.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfpolyx/VCFPolyXTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfpolyx/VCFPolyXTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpolyx** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java  -jar dist/vcfpolyx.jar -R reference.fa input.vcf
(...)
2   1133956 .   A   G   2468.84 .   POLYX=23
2   1133956 .   A   AG  3604.25 .   POLYX=23
2   2981671 .   T   G   47.18   .   POLYX=24
(...)
```

## Cited in:

  * "Multiscale heterogeneity in gastric adenocarcinomaevolution is an obstacle to precision medicine" https://assets.researchsquare.com/files/rs-62554/v1/7883b5d6-a5e6-4d39-8554-e9fef719ac42.pdf
  * Maitena Tellaetxe-Abete, Borja Calvo, Charles Lawrie, Ideafix: a decision tree-based method for the refinement of variants in FFPE DNA sequencing data, NAR Genomics and Bioinformatics, Volume 3, Issue 4, December 2021, lqab092, https://doi.org/10.1093/nargab/lqab092
  * Pol32, an accessory subunit of DNA polymerase delta, plays an essential role in genome stability and pathogenesis of Candida albicans. Shraddheya Kumar Patel & al. Gut Microbes. https://doi.org/10.1080/19490976.2022.2163840 2023.
  * Heczko, L., Hlavac, V., Holy, P. et al. Prognostic potential of whole exome sequencing in the clinical management of metachronous colorectal cancer liver metastases. Cancer Cell Int 23, 295 (2023). https://doi.org/10.1186/s12935-023-03135-x
  * Heczko, L., Liska, V., Vycital, O. et al. Targeted panel sequencing of pharmacogenes and oncodrivers in colorectal cancer patients reveals genes with prognostic significance. Hum Genomics 18, 83 (2024). https://doi.org/10.1186/s40246-024-00644-2


