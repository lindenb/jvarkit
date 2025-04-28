# BedClusterName

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Clusters a BED file into a set of BED files using the 4th column of the bed name.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedclustername  [options] Files

Usage: bedclustername [options] Files
  Options:
    -C, --contig, --chromosome
      group by chromosome.
      Default: false
    -F, --format
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
      Default: BED
      Possible Values: [BED, BED_GZ, INTERVAL_LIST, INTERVAL_LIST_GZ]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -J, --jobs
      number of clusters. (or specify --size or --window-size/--window-shif)
      Default: -1
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    --md5-dir, --sub-dir
      prevent the creation of too many files in the same directory. Create 
      some intermediate directories based on filename's md5.
      Default: false
  * -o, --out
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    -R, --reference
      For sorting and /or writing interval_list,A SAM Sequence dictionary 
      source: it can be a *.dict file, a fasta file indexed with 'picard 
      CreateSequenceDictionary' or 'samtools dict', or any hts file containing 
      a dictionary (VCF, BAM, CRAM, intervals...)
    -S, --size
      number of bases max per bin. (or specify --jobs or 
      --window-size/--window-shif). A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: -1
    --version
      print version and exit

```


## Keywords

 * bed
 * chromosome
 * contig



## Creation Date

2050428

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedclustername/BedClusterName.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedclustername/BedClusterName.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedclustername** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar bedclustername -j 10 -m jeter.mf -o jeter.zip --compress --contig test.bed

```



