# VcfFilterGtf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter VCF on GTF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcffiltergtf  [options] Files

Usage: vcffiltergtf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --biotype
      What biotype 'x'  should be used.
      Default: all
      Possible Values: [all, protein_coding]
    --extend, --extends, -x
      extends each gtf feature by 'x' bases.A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    --gene_id
      Optional file containing a list of gene IDs . Filter on attribute 
      gene_id 
    --gene_name
      Optional file containing a list of gene names . Filter on attribute 
      gene_name 
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -gtf, --gtf
      GTF file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --inverse
      inverse output
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -t, --type
      What should be used in the gtf to keep vcf records
      Default: GENE
      Possible Values: [GENE, TRANSCRIPT, EXON, EXON_BOUDARIES]
    --version
      print version and exit

```


## Keywords

 * vcf
 * gtf



## Creation Date

20230703

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgtf/VcfFilterGtf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgtf/VcfFilterGtf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffiltergtf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Motivation

filter a VCF with a GTF. I'm fed up to filter a  GTF, convert the GTF to bed and use the bed to filter the GTF.


