# BGenToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert bgen to vcf


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bgen2vcf  [options] Files

Usage: bgen2vcf [options] Files
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
    -N, --head
      limit to 'N' variant; negative==show all
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -q, --low-qual
      mark the genotypes for LOWQUAL if all probs are < 'x' and mark the 
      variant if all genotypes are missing or LOW_QUAL. Ignore if <=0
      Default: -1.0
    -G, --no-genotype
      skip genotypes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --targets
      Restrict to comma-separated list of regions. Can be prefixed with "^" to 
      request logical complement
    --targets-file
      Restrict to regions listed in BED 'FILE'. Can be prefixed with "^" to 
      request logical complement.
    -t, --treshold
      Call genotype. Select the best genotype if a probablity is greater that 
      this value. Do not call if the value is <= 0
      Default: -1.0
    --version
      print version and exit

```


## Keywords

 * bgen
 * vcf



## Creation Date

20250508

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bgen/bgen2vcf/BGenToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bgen/bgen2vcf/BGenToVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bgen2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```

```


