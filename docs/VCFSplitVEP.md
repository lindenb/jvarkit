# VCFSplitVEP

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split CSQ vep annotations


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfsplitvep  [options] Files

Usage: vcfsplitvep [options] Files
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
    -o, --out
      Output file. Optional . Default: stdout
    -t, --tag, --tags
      VEP tags name{:type{:aggregate}},name2{:type2{:aggregate2}},etc...  
      where type is a VCF info type Integer,String,Float an aggregate is one 
      of none,min,max,uniq,first
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * vcf
 * vep
 * annotation



## Creation Date

20250517

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsplitvep/VCFSplitVEP.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsplitvep/VCFSplitVEP.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfsplitvep/VCFSplitVEPTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfsplitvep/VCFSplitVEPTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsplitvep** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

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
  *  Critical roles of Dpb3-Dpb4 sub-complex of DNA polymerase epsilon in DNA replication, genome stability, and pathogenesis of Candida albicans. Bhabasha Gyanadeep Utkalaja, Shraddheya Kumar Patel, Satya Ranjan Sahu, Abinash Dutta, Narottam Acharya.  DOI: 10.1128/mbio.01227-24 . mBio. 2024 Aug 29:e0122724.




