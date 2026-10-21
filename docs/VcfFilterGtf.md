# VcfFilterGtf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter/Annotate VCF using a GTF file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfgtf  [options] Files

Usage: vcfgtf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --bed
      only use GTF line overlaping that bed file
    --biotype
      What biotype 'x'  should be used.
      Default: all
      Possible Values: [all, protein_coding]
    --extend, --extends, -x, -X
      extends each gtf feature by 'x' bases. Multiple separated by commas. 
      Empty will be replaced by '0'. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    --gene_id
      Optional comma delimited string or File containing a list of gene IDs . 
      Filter on attribute gene_id. Trailing '.[0-9]+' for version will be 
      removed 
      Default: <empty string>
    --gene_name
      Optional comma delimited string or File containing a list of gene names 
      . Filter on attribute gene_name
      Default: <empty string>
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -gtf, --gtf
      GTF file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-no-overlap, --hide
      Don't show variant that have NO intersection with any of the GTF 
      features 
      Default: false
    --no-biotype
      Do not insert gene biotypes in the output
      Default: false
    --no-gene-id
      Do not insert gene ids in the output
      Default: false
    --no-gene-name
      Do not insert gene names in the output
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
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

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfgtf/VcfFilterGtfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfgtf/VcfFilterGtfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgtf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Motivation

filter/annotate a VCF with a GTF. I'm fed up to filter a  GTF, convert the GTF to bed and use the bed to filter the GTF.


