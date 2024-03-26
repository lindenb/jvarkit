# VcfConcat

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Concatenate VCFs with same sample. See also bcftools concat


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfconcat  [options] Files

Usage: vcfconcat [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --chrom, --contig
      limit to that chromosome
    -G, --drop-genotypes
      Drop genotypes
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --merge
      merge all samples. First Scan all files to get all distinct samples
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -S, --samples
      implies --drop-genotypes
      Default: none
      Possible Values: [none, all, with_alt]
    -T, --tag
      if not empty, add INFO/tag containing the source/path of the variant
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * vcf



## Creation Date

20131230

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfConcat.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfconcat/VcfConcat.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfconcat** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

concat numerous VCF, without checking files at start, doesn't sort.

## Input

input is a set of VCF files or a file with the suffix .list containing the path to the vcfs

## Example

### From stdin

```bash
$ find ./ -name "*.vcf" | grep Sample1 | java -jar dist/jvarkit vcfconcat > out.vcf
```

### From files

```bash
$ java -jar dist/jvarkit vcfconcat Sample1.samtools.vcf Sample1.gatk.vcf > out.vcf
$ java -jar dist/jvarkit vcfconcat paths.list > out.vcf
```


## See also

bcftools concat


