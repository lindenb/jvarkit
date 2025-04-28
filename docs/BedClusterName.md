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
    --merge-distance, --merge
      if greater than -1 Merge overlapping record for the same name within a 
      distance of 'x'. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
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
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | awk -F '\t' '$3=="exon"' | \
	java -jar dist/jvarkit.jar gtf2bed -c 'gene_name' |\
	java -jar dist/jvarkit.jar bedclustername -o jeter.zip  --size 100 --merge 1

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     7332  2025-04-28 14:37   cluster.000000001.bed
     1740  2025-04-28 14:37   cluster.000000002.bed
     1456  2025-04-28 14:37   cluster.000000003.bed
---------                     -------
    10528                     3 files


$ unzip -p jeter.zip cluster.000000002.bed | head
22	41697718	41697776	ZC3H7B
22	41716664	41716717	ZC3H7B
22	41721567	41721601	ZC3H7B
22	41721724	41721922	ZC3H7B
22	41723209	41723368	ZC3H7B
22	41726026	41726107	ZC3H7B
22	41728174	41732847	ZC3H7B
22	41734316	41734359	ZC3H7B
22	41735004	41735195	ZC3H7B
22	41735819	41736141	ZC3H7B

```



