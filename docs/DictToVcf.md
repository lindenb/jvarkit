# DictToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert a SAM dictionary from vcf,sam,bam,dict, etc.. to vcf.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar dict2vcf  [options] Files

Usage: dict2vcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --extra-lines
      File containing extra header lines
    --format, --formats
      Add those predefined FORMAT fields
      Default: GT,GQ,DP,PL,AD,FT,PS,PQ
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --info, --infos
      Add those predefined INFO fields
      Default: END,DB,DP,SB,AF,AC,AN,MQ0,MQ,SOMATIC
    -o, --out
      Output file. Optional . Default: stdout
    --samples
      Add those samples
      Default: []
    --samples-file
      Read samples from file. If the file has a BAM/CRAM/BCF/VCF extension, 
      samples are extracted from their metadata.
    --version
      print version and exit

```


## Keywords

 * dict
 * bed
 * sam
 * bam
 * vcf



## Creation Date

20251126

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2vcf/DictToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dict2vcf/DictToVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/dict2vcf/DictToVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/dict2vcf/DictToVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **dict2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

convert Dict to VCF when you want to generate an empty vcf file

## Example

```
java -jar dist/jvarkit.jar dict2vcf src/test/resources/human_b37.dict  --samples S1,S3
```



