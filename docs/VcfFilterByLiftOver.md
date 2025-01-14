# VcfFilterByLiftOver

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Add FILTER(s) to a variant when it is known to map elsewhere after liftover.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcffilterbyliftover  [options] Files

Usage: vcffilterbyliftover [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
  * -f, --chain
      LiftOver chain file. Can be a local chain file, a URL 'https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToCriGri1.over.chain.gz', 
      or a chain identifier like 'hg19ToHg38'.
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --lifted-filtered-vcf
      Another VCF in the destination reference. Input VCF will be filtered if 
      this vcf is FILTERed at the same lifted position.
    -m, --minmatch
      lift over min-match.
      Default: 1.0
    --no-validation
      Disable dictionary validation
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -D
      final distance. We look for weird distance between the current and the 
      previous variant on the same contig.Two variants initially distance < d 
      should have a distance <D after lift over.
      Default: 1500
    -d
      initial distance. See option -D.
      Default: 1000

```


## Keywords

 * vcf
 * liftover



## Creation Date

20190418

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfFilterByLiftOver.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/VcfFilterByLiftOver.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffilterbyliftover** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


