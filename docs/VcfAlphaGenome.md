# VcfAlphaGenome

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate VCF with Tabix index files from Google AlphaGenome


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfalphagenome  [options] Files

Usage: vcfalphagenome [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --buffer-size
      Buffer size / processing window size. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 1000
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --impact
      Path to tabix indexed impact_score_snvs.tsv.gz
    -na, --na
      value  used for 'allele-not-found'.
      Default: -999.0
    -o, --out
      Output file. Optional . Default: stdout
    --splicing, --combined
      Path to tabix indexed combined_splicing_snvs.tsv.gz
    --tag
      INFO tag prefix.
      Default: AG_
    --version
      print version and exit

```


## Keywords

 * vcf
 * prediction
 * google
 * alphagenome



## Creation Date

20260915

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfalphagenome/VcfAlphaGenome.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfalphagenome/VcfAlphaGenome.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfalphagenome** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example
 
 ```
 $ bcftools view in.vcf |\
   java -jar dist/jvarkit.jar vcfalphagenome \
 	--impact   /path/to/google/AlphaGenome/alphagenome_variant_impact_score_snvs.GRCh38.tsv.gz \\
 	--splicing /path/to/google/AlphaGenome/combined_alphagenome_splicing_snvs.GRCh38.tsv.gz |\
 	grep -v "##" | head -n 4

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    26169987        .       G       A       .       .       AG_PHRED=22.68624;AG_RAW=0.7722;AG_SPLICING=0.08169
chr1    42759328        .       G       A       .       .       AG_PHRED=18.73157;AG_RAW=0.4941;AG_SPLICING=0.06523
chr1    43558671        .       C       T       .       .       AG_PHRED=14.63943;AG_RAW=0.307;AG_SPLICING=0.05377
```
 

