# ShiftVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

shit all coordinates of a VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar shiftvcf  [options] Files

Usage: shiftvcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --bed
      Bed Path. Extract Variants overlaping this BED. Or use -R2.
    -R2, --destination-reference
      Original fasta reference. We shift the VCF back to this reference. 
      Required without --bed.A SAM Sequence dictionary source: it can be a 
      *.dict file, a fasta file indexed with 'picard CreateSequenceDictionary' 
      or 'samtools dict', or any hts file containing a dictionary (VCF, BAM, 
      CRAM, intervals...)
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * bed



## Creation Date

20241002

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/shiftvcf/ShiftVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/shiftvcf/ShiftVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **shiftvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




shift coordinates of a VCF.


## Usage 1

shift a vcf using a using ROI.bed

```
bcftools view --regions-file ROI.bed source.vcf.gz |\
	java -jar dist/jvarkit.jar shiftvcf --bed ROI.bed
```


## Usage 2 

Shift the VCF back to original coordinates


```
java -jar dist/jvarkit.jar shiftvcf --R2 source.vcf.gz shifted.vcf

````





