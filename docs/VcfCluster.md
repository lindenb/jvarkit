# VcfCluster

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfcluster  [options] Files

Usage: vcfcluster [options] Files
  Options:
    -f, --force
      force writing
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -L, --length, -n
      max-length per cluster. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 1000000
    --no-group-bnd
      disable grouping variants with INFO/SVTYPE=BND  having same 
      ID=INFO/MATEID 
      Default: false
    --no-index
      do no write tbi index
      Default: false
  * -o, --output
      output directory
    -p, --prefix
      file prefix
      Default: split.
    -seed, --seed
      random see: -1= current-time. Random is used to add a variant to a 
      random available batch
      Default: -1
    --version
      print version and exit

```


## Keywords

 * vcf



## Creation Date

20260108

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcluster/VcfCluster.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcluster/VcfCluster.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcluster/VcfClusterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcluster/VcfClusterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcluster** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

split CNV/SV vcf by clusters of size='x' for later parallelization (e.g: duphold, annotation, etc...)
All INFO/SVTYPE=BND are grouped in the same cluster using ID and INFO/MATEID

## Example

```bash
java -jar dist/jvarkit.jar vcfcluster -o TMP manta.vcf.gz

```


